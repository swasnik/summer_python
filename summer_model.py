
import numpy


class EpiModel:
    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 report=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 default_starting_compartment="", equilibrium_stopping_tolerance=None, output_connections=(),
                 tracked_quantities=()):

        # attributes with specific format that are independent of user inputs
        self.compartment_values = {}
        self.transition_flows, self.death_flows = [[]] * 2

        # features that should not be changed
        self.available_birth_approaches = ["add_crude_birth_rate", "replace_deaths", "no_births"]

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
            infectious_compartment, birth_approach, report, reporting_sigfigs, entry_compartment,
            starting_population, default_starting_compartment, equilibrium_stopping_tolerance, output_connections,
            tracked_quantities)

        # stop ide complaining about attributes being defined outside __init__, even though they aren't
        self.times, self.compartment_types, self.initial_conditions, self.parameters, self.requested_flows, \
            self.initial_conditions_to_total, self.infectious_compartment, self.birth_approach, self.report, \
            self.reporting_sigfigs, self.entry_compartment, self.starting_population, \
            self.default_starting_compartment, self.default_starting_population, self.equilibrium_stopping_tolerance, \
            self.output_connections, self.tracked_quantities, self.unstratified_flows = [None] * 18

        # convert input arguments to model attributes
        for attribute in ["times", "compartment_types", "initial_conditions", "parameters",
                          "initial_conditions_to_total", "infectious_compartment", "birth_approach", "report",
                          "reporting_sigfigs", "entry_compartment", "starting_population",
                          "default_starting_compartment", "infectious_compartment", "equilibrium_stopping_tolerance",
                          "output_connections", "tracked_quantities"]:
            setattr(self, attribute, eval(attribute))

        # set initial conditions and implement flows
        self.set_initial_conditions(initial_conditions_to_total)

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.add_default_quantities()

    def output_to_user(self, comment):
        """
        short function to save the if statement in every call to output some information
        """
        if self.report:
            print(comment)

    def check_and_report_attributes(
            self, times, compartment_types, initial_conditions, parameters, requested_flows,
            initial_conditions_to_total, infectious_compartment, birth_approach, report, reporting_sigfigs,
            entry_compartment, starting_population, default_starting_compartment, equilibrium_stopping_tolerance,
            output_connections, tracked_quantities):
        """
        check all input data have been requested correctly
        """

        # check that variables are of the expected type
        for expected_numeric_variable in ["reporting_sigfigs", "starting_population"]:
            if not isinstance(eval(expected_numeric_variable), int):
                raise TypeError("expected integer for %s" % expected_numeric_variable)
        for expected_list in ["times", "compartment_types", "requested_flows"]:
            if not isinstance(eval(expected_list), list):
                raise TypeError("expected list for %s" % expected_list)
        for expected_string in \
                ["infectious_compartment", "birth_approach", "entry_compartment", "default_starting_compartment"]:
            if not isinstance(eval(expected_string), str):
                raise TypeError("expected string for %s" % expected_string)
        for expected_boolean in ["initial_conditions_to_total", "report"]:
            if not isinstance(eval(expected_boolean), bool):
                raise TypeError("expected boolean for %s" % expected_boolean)
        for expected_dict in ["output_connections", "tracked_quantities"]:
            if not isinstance(eval(expected_dict), tuple):
                raise TypeError("expected dictionary for %s" % expected_dict)

        # check some specific requirements
        if infectious_compartment not in compartment_types:
            ValueError("infectious compartment name is not one of the listed compartment types")
        if birth_approach not in self.available_birth_approaches:
            ValueError("requested birth approach unavailable")
        if sorted(times) != times:
            self.output_to_user("requested integration times are not sorted, now sorting")
            self.times = sorted(self.times)

        # report on characteristics of inputs
        if report:
            print("integrating from time %s to %s"
                  % (round(times[0], reporting_sigfigs), round(times[-1], reporting_sigfigs)))
            print("unstratified requested initial conditions are:")
            for compartment in initial_conditions:
                print("\t%s: %s" % (compartment, initial_conditions[compartment]))
            print("infectious compartment is called '%s'" % infectious_compartment)
            print("birth approach is %s" % birth_approach)

    def set_initial_conditions(self, initial_conditions_to_total):
        """
        set starting compartment values
        """
        # set starting values of unstratified compartments to requested value, or zero if no value requested
        for compartment in self.compartment_types:
            if compartment in self.initial_conditions:
                self.compartment_values[compartment] = self.initial_conditions[compartment]
            else:
                self.output_to_user("no starting value requested for %s so set to zero" % compartment)
                self.compartment_values[compartment] = 0

        # sum to a total value if requested
        if initial_conditions_to_total:
            self.sum_initial_compartments_to_total()

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain value
        """
        compartment = self.find_remainder_compartment()
        if sum(self.compartment_values.values()) > self.starting_population:
            raise ValueError("total of requested compartment values is greater than the requested starting population")
        remaining_population = self.starting_population - sum(self.compartment_values.values())
        self.output_to_user("requested that total population sum to %s" % self.starting_population)
        self.output_to_user("remaining population of %s allocated to %s compartment"
                            % (remaining_population, compartment))
        self.compartment_values[compartment] = remaining_population

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total
        """
        if len(self.default_starting_compartment) > 0 and \
                self.default_starting_compartment not in self.compartment_values:
            raise ValueError("starting compartment to populate with initial values not found in available compartments")
        elif len(self.default_starting_compartment) > 0:
            return self.default_starting_compartment
        else:
            self.output_to_user("no default starting compartment requested for unallocated population, " +
                                "so will be allocated to entry compartment %s" % self.entry_compartment)
            return self.entry_compartment

    def implement_flows(self, requested_flows):
        """
        add all flows to create data frames from input lists
        """
        for flow in requested_flows:
            if flow["parameter"] not in self.parameters:
                raise ValueError("flow parameter not found in parameter list")
            if flow["from"] not in self.compartment_types:
                raise ValueError("from compartment name not found in compartment types")
            if "to" in flow and flow["to"] not in self.compartment_types:
                raise ValueError("to compartment name not found in compartment types")

            if flow["type"] == "compartment_death":
                self.add_death_flow(flow)
            else:
                self.add_transition_flow(flow)

        # retain a copy of the original flows for the purposes of graphing, etc.
        self.unstratified_flows = self.transition_flows

    def add_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed
        """

        pass

    def add_transition_flow(self, flow):
        """
        simply add a flow to the data frame storing the flows
        """
        flow["implement"] = 0
        self.transition_flows.append(flow)

    def add_death_flow(self, flow):
        """
        similarly for compartment-specific death flows
        """
        flow["implement"] = 0
        self.death_flows.append(flow)


if __name__ == "__main__":
    sir_model = EpiModel(numpy.linspace(0, 60 / 365, 60).tolist(),
                         ["susceptible", "infectious", "recovered"],
                         {"infectious": 0.001},
                         {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
                         [{"type": "standard_flows", "parameter": "recovery", "from": "infectious", "to": "recovered"},
                          {"type": "infection_density", "parameter": "beta", "from": "susceptible", "to": "infectious"},
                          {"type": "compartment_death", "parameter": "infect_death", "from": "infectious"}],
                         report=False)
    print(sir_model.death_flows)
