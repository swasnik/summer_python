
import numpy


class EpiModel:
    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 report_progress=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 default_starting_compartment="", equilibrium_stopping_tolerance=None, output_connections=(),
                 tracked_quantities=()):

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
            infectious_compartment, birth_approach, report_progress, reporting_sigfigs, entry_compartment,
            starting_population, default_starting_compartment, equilibrium_stopping_tolerance, output_connections,
            tracked_quantities)

        # convert input arguments to model attributes
        for attribute in ["times", "compartment_types", "initial_conditions", "parameters", "infectious_compartment",
                          "birth_approach", "report_progress", "reporting_sigfigs", "entry_compartment",
                          "starting_population", "default_starting_compartment", "infectious_compartment",
                          "equilibrium_stopping_tolerance", "output_connections", "tracked_quantities"]:
            setattr(self, attribute, eval(attribute))

        # set initial conditions and implement flows
        self.set_initial_conditions(initial_conditions_to_total)

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.add_default_quantities()

    def check_and_report_attributes(
            self, times, compartment_types, initial_conditions, parameters, requested_flows,
            initial_conditions_to_total, infectious_compartment, birth_approach, report_progress, reporting_sigfigs,
            entry_compartment, starting_population, default_starting_compartment, equilibrium_stopping_tolerance,
            output_connections, tracked_quantities):
        """
        check all input data are in the correct form
        """

        pass

    def set_initial_conditions(self, initial_conditions_to_total):
        """
        set starting compartment values
        """

        pass

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain value
        """

        pass

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total
        """

        pass

    def implement_flows(self, requested_flows):
        """
        add all flows to create data frames from input lists
        """

        pass

    def add_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed
        """

        pass


if __name__ == "__main__":
    sir_model = EpiModel(numpy.linspace(0, 60 / 365, 60),
                         ["susceptible", "infectious", "recovered"],
                         {"infectious": 0.001},
                         {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
                         [{"type": "standard_flows", "parameter": "recovery", "from": "infectious", "to": "recovered"},
                          {"type": "infection_density", "parameter": "beta", "from": "susceptible", "to": "infectious"},
                          {"type": "compartment_death", "parameter": "infect_death", "from": "infectious"}])
