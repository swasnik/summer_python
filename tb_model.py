
import summer_model
import numpy
import matplotlib.pyplot


def add_vtp_latency_parameters(parameters_, age_group, change_time_unit=365.25):
    """
    function to add the latency parameters estimated by Ragonnet et al from our paper in Epidemics to the existing
    parameter dictionary
    """
    vtp_latency_parameters = \
        {"all_ages":
            {"early_progression": 1.1e-3,
             "stabilisation": 1.0e-2,
             "late_progression": 5.5e-6},
         0: {"early_progression": 6.6e-3,
             "stabilisation": 1.2e-2,
             "late_progression": 1.9e-11},
         5: {"early_progression": 2.7e-3,
             "stabilisation": 1.2e-2,
             "late_progression": 6.4e-6},
         15: {"early_progression": 2.7e-4,
              "stabilisation": 5.4e-3,
              "late_progression": 3.3e-6}}
    latency_parameters = \
        {key: value * change_time_unit for key, value in vtp_latency_parameters[age_group].items()}
    parameters_.update(latency_parameters)
    return parameters_


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
        {"beta": 4.0,
         "recovery": case_fatality_rate / untreated_disease_duration,
         "infect_death": (1.0 - case_fatality_rate) / untreated_disease_duration,
         "universal_death_rate": 1.0 / 50.0}
    parameters = add_vtp_latency_parameters(parameters, "all_ages")
    times = numpy.linspace(0.0, 200.0, 200).tolist()
    flows = [{"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "early_latent"},
             {"type": "infection_frequency", "parameter": "beta", "origin": "recovered", "to": "early_latent"},
             {"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
             {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}]
    flows = add_standard_latency_flows(flows)

    # instantiate and run
    tb_model = summer_model.StratifiedModel(
        times, ["susceptible", "early_latent", "late_latent", "infectious", "recovered"], {"infectious": 1e-3},
        parameters, flows, birth_approach="replace_deaths", report=True)
    tb_model.run_model()

    # get outputs
    matplotlib.pyplot.plot(times, tb_model.outputs[:, tb_model.compartment_names.index("infectious")] * 1e5)
    # matplotlib.pyplot.xlim((1e3, 2e3))
    # matplotlib.pyplot.ylim((0.0, 100.0))
    matplotlib.pyplot.show()
