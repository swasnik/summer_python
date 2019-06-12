
import summer_model
import numpy
import matplotlib.pyplot


def add_vtp_latency_parameters(parameters_, change_time_unit=365.25):
    """
    function to add the latency parameters estimated by Ragonnet et al from our paper in Epidemics to the existing
    parameter dictionary
    """
    vtp_latency_parameters = \
        {"early_progression": 1.1e-3,
         "stabilisation": 1.0e-2,
         "late_progression": 5.5e-6}
    parameters_.update({key: value * change_time_unit for key, value in vtp_latency_parameters.items()})
    return parameters_


def get_age_specific_latency_parameters(parameter, change_time_unit=365.25):
    age_stratified_parameters = \
        {"early_progression":
             {"0": 6.6e-3,
              "5": 2.7e-3,
              "15": 2.7e-4},
         "stabilisation":
             {"0": 1.2e-2,
              "5": 1.2e-2,
              "15": 5.4e-3},
         "late_progression":
             {"0": 1.9e-11,
              "5": 6.4e-6,
              "15": 3.3e-6}}
    result = {key: value * change_time_unit for key, value in age_stratified_parameters[parameter].items()}
    return result


def add_standard_latency_flows(flows_):
    """
    adds our standard latency flows to the list of flows to be implemented in the model
    """
    flows_ += [
        {"type": "standard_flows", "parameter": "early_progression", "origin": "early_latent", "to": "infectious"},
        {"type": "standard_flows", "parameter": "stabilisation", "origin": "early_latent", "to": "late_latent"},
        {"type": "standard_flows", "parameter": "late_progression", "origin": "late_latent", "to": "infectious"}]
    return flows_


if __name__ == "__main__":

    # set basic parameters, flows and times, except for latency flows and parameters, then functionally add latency
    case_fatality_rate = 0.4
    untreated_disease_duration = 3.0
    parameters = \
        {"beta": 50.0,
         "recovery": case_fatality_rate / untreated_disease_duration,
         "infect_death": (1.0 - case_fatality_rate) / untreated_disease_duration,
         "universal_death_rate": 1.0 / 50.0}
    parameters = add_vtp_latency_parameters(parameters)
    times = numpy.linspace(0.0, 200.0, 200).tolist()
    flows = [{"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "early_latent"},
             {"type": "infection_frequency", "parameter": "beta", "origin": "recovered", "to": "early_latent"},
             {"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
             {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}]
    flows = add_standard_latency_flows(flows)

    # instantiate and run
    tb_model = summer_model.StratifiedModel(
        times, ["susceptible", "early_latent", "late_latent", "infectious", "recovered"], {"infectious": 1e-3},
        parameters, flows, birth_approach="replace_deaths")
    tb_model.stratify("age", [5, 15], [],
                      adjustment_requests=
                      {"early_progression": {"adjustments": get_age_specific_latency_parameters("early_progression")},
                       "stabilisation": {"adjustments": get_age_specific_latency_parameters("stabilisation")},
                       "late_progression": {"adjustments": get_age_specific_latency_parameters("late_progression")}},
                      report=True)
    tb_model.run_model()

    # get outputs
    matplotlib.pyplot.plot(times, tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_15")] * 1e5)
    # matplotlib.pyplot.xlim((1e3, 2e3))
    # matplotlib.pyplot.ylim((0.0, 100.0))
    matplotlib.pyplot.show()
