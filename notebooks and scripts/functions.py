import re
import pandas as pd
import numpy as np
import cobra
import matplotlib.pyplot as plt
from cobra import Reaction, Metabolite


def retag_model_rxns_and_mets(model, tag):
    """
    create copy of the model with all metabolites and reactions with a tag at the end of the id
    """
    
    temp_model = model.copy()
    
    for met in temp_model.metabolites:
        met.id = met.id + tag
        met.compartment = met.compartment + tag
    for rxn in temp_model.reactions:
        rxn.id = rxn.id + tag
        
    return temp_model


###################### functions to set ratio constraints ############################################################

# [] make more reusable by changing "_00" and "_12" for user-defined day and night tags

def set_output_day_night_ratio(model, ratio, output_string = "Marchantia_biomass"):
    """
    Sets the ratio between the daytime output reaction and night-time output reaction.
    """
    try:
        model.remove_cons_vars(model.constraints.output_light_dark_const)
        print("output_light_dark_const constraint removed")
    except:
        #print("output_light_dark_const not found")
        pass
        
    output_light_dark_const = model.problem.Constraint(
                            (model.reactions.get_by_id(output_string + "_light").flux_expression) - 
                            ratio*(model.reactions.get_by_id(output_string + "_dark").flux_expression), 
                            lb = 0,
                            ub= 0,
                            name="output_light_dark_const")
    model.add_cons_vars(output_light_dark_const)
                          
    return ratio, ""


def set_Rubisco_C_to_O_ratio(model, ratio):
    """
    sets the ratio between Rubisco carboxylation and oxygenation
    0.1 => 1:10 C:O
    10 => 10:1 C:O
    """
    try:
        model.remove_cons_vars(model.constraints.rubisco_c_to_o_const00)
        print("rubisco_c_to_o_const00 constraint removed")

    except:
        #print("rubisco_c_to_o_const00 not found")
        pass

    rubisco_c_to_o_const00 = model.problem.Constraint(
                            (model.reactions.RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_00.flux_expression) - 
                            ratio*(model.reactions.RXN_961_p_00.flux_expression), 
                            lb = 0,
                            ub= 0,
                            name="rubisco_c_to_o_const00")
    model.add_cons_vars(rubisco_c_to_o_const00)
    return ratio, ""


def set_NGAM_ATP_NADPH_ratio(model, ratio):
    """
    Set constraint between ATPase reaction and sum of NADPH oxidases reactions
    """
    try:
        model.remove_cons_vars(model.constraints.ATPase_to_NADPHox_const00)
        print("ATPase_to_NADPHox_const00 constraint removed")
    except:
        #print("ATPase_to_NADPHox_const00 not found")
        pass
    try:
        model.remove_cons_vars(model.constraints.ATPase_to_NADPHox_const12)
        print("ATPase_to_NADPHox_const12 constraint removed")
    except:
        pass
    
    for time in ["00", "12"]:        

        ATPase_NADPHox_const = model.problem.Constraint(
            (model.reactions.get_by_id(f"ATPase_tx_{time}").flux_expression) - 
            ratio*3*model.reactions.get_by_id(f"NADPHoxc_tx_{time}").flux_expression, 
            lb=0,
            ub=0,
            name=f"ATPase_to_NADPHox_const{time}"
        )
        model.add_cons_vars(ATPase_NADPHox_const)

        for compartment in ["m", "p"]:
            NADPHox_compartment_const = model.problem.Constraint(
                (model.reactions.get_by_id(f"NADPHoxc_tx_{time}").flux_expression) - 
                model.reactions.get_by_id(f"NADPHox{compartment}_tx_{time}").flux_expression, 
                lb=0,
                ub=0,
                name=f"NADPHox_{compartment}_const{time}"
            )
            model.add_cons_vars(NADPHox_compartment_const)

    return ratio, ""


def set_NGAM_day_night_ratio(model, ratio):
    """
    set up constraint for ratio between ATPase day and ATPase_night activity to set NGAM day:night constraint
    3 => 3:1 day: night
    """
    
    try:
        model.remove_cons_vars(model.constraints.ATPase_light_dark)
        print("ATPase_light_dark constraint removed")
    except:
        #print("ATPase_light_dark not found")
        pass
    
    ATPase_light_dark = model.problem.Constraint(
                            (model.reactions.get_by_id("ATPase_tx_00").flux_expression) - 
                            ratio*(model.reactions.get_by_id("ATPase_tx_12").flux_expression), 
                            lb = 0,
                            ub= 0,
                            name="ATPase_light_dark")
    model.add_cons_vars(ATPase_light_dark)
    return ratio, ""


def set_nitrate_uptake_day_night_ratio(model, ratio):
    """
    set constraint between nitrate uptake during the day and nitrate uptake during the night
    """
    try:
        model.remove_cons_vars(model.constraints.nitrate_light_dark_const)
        print("nitrate_light_dark_const constraint removed")
    except:
        #print("nitrate_light_dark_const not found")
        pass

    nitrate_light_dark_const = model.problem.Constraint(
                            (model.reactions.Nitrate_tx_00.flux_expression) - 
                            ratio*(model.reactions.Nitrate_tx_12.flux_expression), 
                            lb = 0,
                            ub= 0,
                            name= "nitrate_light_dark_const")
    model.add_cons_vars(nitrate_light_dark_const)
    
    return ratio, ""

def set_photorespiratory_loss_ratio(model,reassim_resp_CO2):
    """
    set constraint between CO2_m_00 (mitochondrial CO2 loss reaction) and CO2_mc_00 (mitochondrial CO2 export to cytosol)
    """
    try:
        model.remove_cons_vars(model.constraints.CO2_resp_loss_const)
        print("CO2_resp_loss_const constraint removed")
    except:
        print("CO2_resp_loss_const not found")
        pass
    

    CO2_resp_loss_const = model.problem.Constraint(
                            (model.reactions.CO2_m_00.flux_expression) - 
                            (1-reassim_resp_CO2)*(model.reactions.CO2_mc_00.flux_expression), 
                            lb = 0,
                            ub= 0,
                            name= "CO2_resp_loss_const")
    model.add_cons_vars(CO2_resp_loss_const)
    
    return CO2_resp_loss_const, ""

#############################################################################################################

def model_results_overview(model, solut, objective_reaction, rounding_nr=5, linker_tag="linker", light_tag="_00", dark_tag="_12", linker_tag1="_00_to_12", linker_tag2="_12_to_00"):
    linker_flux_dict = {}
    for rxn in model.reactions:
        if linker_tag in rxn.id:
            linker_id = rxn.id[:-9]
            if linker_id not in linker_flux_dict:
                linker_flux_dict[linker_id] = [0, 0]
            if rxn.id.endswith(linker_tag1):
                linker_flux_dict[linker_id][0] = round(rxn.flux, rounding_nr)
            elif rxn.id.endswith(linker_tag2):
                linker_flux_dict[linker_id][1] = round(rxn.flux, rounding_nr)

    sum_of_fluxes = 0
    for flux in solut:
        sum_of_fluxes += abs(flux)
    


    print(f"Photon uptake = {round(model.reactions.get_by_id('Photon_tx'+light_tag).flux, 4)}   % of allowed Photon uptake = {round((model.reactions.get_by_id('Photon_tx'+light_tag).flux/model.reactions.get_by_id('Photon_tx'+light_tag).upper_bound)*100, 2)}")
    print(f"Output rate ({objective_reaction}) {round(model.reactions.get_by_id(objective_reaction).flux, 4)}\t\tsum of fluxes: {round(sum_of_fluxes, 4)}")
    print(f"gas exchange = Day: {round(model.reactions.get_by_id('CO2_tx'+light_tag).flux, 5)} Night: {round(model.reactions.get_by_id('CO2_tx'+dark_tag).flux, 5)}")
    if model.reactions.get_by_id('CO2_tx'+dark_tag).flux <= 0:
        CCE = calculate_CCE(model)
        print(f"CCE: {CCE}")
    else:
        print("CO2_tx at night is taking up CO2.")
    print(f"ATPase: {round(model.reactions.get_by_id('ATPase_tx'+light_tag).flux, 4)} {round(model.reactions.get_by_id('ATPase_tx'+dark_tag).flux, 4)}")
    print(f"Rubisco Carbox./Oxygen. = {round(model.reactions.get_by_id('RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p'+light_tag).flux, 4)}/{round(model.reactions.get_by_id('RXN_961_p'+light_tag).flux, 4)}")
    # print(f"Photorespiratory CO2 losses = {round(model.reactions.get_by_id('CO2_m'+light_tag).flux, 4)} / mitochondrial CO2 export = {round(model.reactions.get_by_id('CO2_mc'+light_tag).flux, 4)}")
    print("\nLinker fluxes \t\t Day \t Night")
    for linker_id, fluxes in linker_flux_dict.items():
        day_flux = "\033[91m" + str(fluxes[0]) + "\033[00m" if round(fluxes[0], rounding_nr) == round(model.reactions.get_by_id(linker_id+linker_tag1).upper_bound, rounding_nr) else str(fluxes[0])
        night_flux = "\033[91m" + str(fluxes[1]) + "\033[00m" if round(fluxes[1], rounding_nr) == round(model.reactions.get_by_id(linker_id+linker_tag2).upper_bound, rounding_nr) else str(fluxes[1])
        print(f"{linker_id:<20}{day_flux:>20}{night_flux:>20}\n")
    calc_carbon_flux(model, output_string=objective_reaction[:-6])


def calculate_CCE(model):

    CO2_day = model.reactions.get_by_id("CO2_tx_00").flux
    CO2_night = model.reactions.get_by_id("CO2_tx_12").flux
    CO2_resp = model.reactions.get_by_id("CO2_m_00").flux if "CO2_m_00" in model.reactions else 0

    
    if CO2_night <= 0:
        try:
            CCE = 1 - ((abs(CO2_night)+abs(CO2_resp))/CO2_day)
        except ZeroDivisionError:
            CCE = 0
            print("CO2_tx_00 = 0")
    else:
        CCE = 1
        
    return CCE

##

def calc_carbon_flux(model, output_string="Marchantia_biomass", night_tag = "_12"):
    """
    Calculates various carbon fluxes and returns the results as a pandas dataframe.
    """
    # Calculate linker carbon fluxes
    day_to_night_carbons = sum(rxn.flux * rxn.reactants[0].elements.get("C", 0) for rxn in model.reactions if "linker_00_to_12" in rxn.id)
    night_to_day_carbons = sum(rxn.flux * rxn.reactants[0].elements.get("C", 0) for rxn in model.reactions if "linker_12_to_00" in rxn.id)
    
    # Print CO2 exchange and respiratory losses
    print("CO2 daytime exchange:", model.reactions.CO2_tx_00.flux)
    night_CO2 = model.reactions.get_by_id("CO2_tx"+night_tag).flux
    print("CO2 nightime exchange:", night_CO2)
    print("CO2 respiratory losses:", model.reactions.CO2_m_00.flux if "CO2_m_00" in model.reactions else "0 (no CO2_m_00 found)")
    
    # Print carbon fluxes
    print("Carbon day to night:", day_to_night_carbons)
    print("Carbon night to day:", night_to_day_carbons)
    print("CO2 exchange night:day:", night_to_day_carbons / day_to_night_carbons)
    
    # Print carbon biomass fluxes
    carbon_biomass_dict = {}
    for phase in ("_00", "_12"):
        carbon_stoich = sum(met.elements.get("C", 0) * model.reactions.get_by_id(output_string + phase).get_coefficient(met.id) for met in model.reactions.get_by_id(output_string + phase).reactants if "C" in met.elements)
        carbon_biomass_dict[model.reactions.get_by_id(output_string + phase).id] = carbon_stoich * model.reactions.get_by_id(output_string + phase).flux
    for item in carbon_biomass_dict:
        print(item, "carbon flux:", carbon_biomass_dict[item])
    print("total carbon biomass flux:", sum(carbon_biomass_dict.values()))
    
    # Print CCE and nocturnal CO2 refixation coefficient
    CCE = calculate_CCE(model)
    print("total CCE:", CCE)
    nCO2RC = calculate_noct_CO2_refixation_coefficient(model, verbose=False, threshold=0.00001, show_noflux=False)
    print("nCO2RC:", nCO2RC)
    
    # Print night CCE
    cce_night = abs(carbon_biomass_dict[model.reactions.get_by_id(output_string + "_12").id]) / day_to_night_carbons
    print("night CCE1 (carbon in to biomass):", cce_night)
    cce_night2 =  1 + (model.reactions.CO2_tx_12.flux / day_to_night_carbons)
    print("night CCE2 (1 - (carbon lost / carbon in)):", cce_night2)

    carbon_flux_dict = {"CO2 daytime exchange": model.reactions.CO2_tx_00.flux,
                        "CO2 nightime exchange": model.reactions.CO2_tx_12.flux,
                        "Carbon day to night": day_to_night_carbons,
                        "Carbon night to day": night_to_day_carbons,
                        "CO2 exchange night:day": night_to_day_carbons / day_to_night_carbons,
                        "day carbon biomass flux": abs(carbon_biomass_dict[output_string + "_00"]),
                        "night carbon biomass flux": abs(carbon_biomass_dict[output_string + "_12"]),
                        "total Carbon biomass flux": abs(sum(carbon_biomass_dict.values())),
                        "total CCE": CCE,
                        "nCO2RC": nCO2RC,
                        "% nCO2R": nCO2RC * 100,
                        "night CCE1": cce_night,
                        "night CCE2: ": cce_night2}

    # create a pandas dataframe from the dictionary
    carbon_flux_df = pd.DataFrame.from_dict(carbon_flux_dict, orient='index', columns=['value'])

    return carbon_flux_df
    
def calc_carbon_flux_single_phase(model, output_string="Marchantia_biomass"):
    """
    Calculates various carbon fluxes and returns the results as a pandas dataframe.
    """
    # Calculate linker carbon fluxes
    day_to_night_carbons = sum(rxn.flux * rxn.reactants[0].elements.get("C", 0) for rxn in model.reactions if ("accumulation" in rxn.id and rxn.flux > 0))
    night_to_day_carbons = sum(rxn.flux * rxn.reactants[0].elements.get("C", 0) for rxn in model.reactions if ("accumulation" in rxn.id and rxn.flux < 0))
    
    # Print CO2 exchange and respiratory losses
    print("CO2 exchange:", model.reactions.CO2_tx.flux)
    print("CO2 respiratory losses:", model.reactions.CO2_m.flux if "CO2_m" in model.reactions else "0 (no CO2_m found)")
    
    # Print carbon fluxes
    print("Carbon day to night:", day_to_night_carbons)
    print("Carbon night to day:", night_to_day_carbons)
    # print("CO2 exchange night:day:", night_to_day_carbons / day_to_night_carbons)
    
    # # Print carbon biomass fluxes
    # carbon_biomass_dict = {}
    # for phase in ("_00", "_12"):
    #     carbon_stoich = sum(met.elements.get("C", 0) * model.reactions.get_by_id(output_string + phase).get_coefficient(met.id) for met in model.reactions.get_by_id(output_string + phase).reactants if "C" in met.elements)
    #     carbon_biomass_dict[model.reactions.get_by_id(output_string + phase).id] = carbon_stoich * model.reactions.get_by_id(output_string + phase).flux
    # for item in carbon_biomass_dict:
    #     print(item, ":", carbon_biomass_dict[item])
    # print("total Carbon biomass flux:", sum(carbon_biomass_dict.values()))
    
    # # Print CCE and nocturnal CO2 refixation coefficient
    # CCE = calculate_CCE(model)
    # print("total CCE:", CCE)
    # nCO2RC = calculate_noct_CO2_refixation_coefficient(model, verbose=False, threshold=0.00001, show_noflux=False)
    # print("nCO2RC:", nCO2RC)
    
    # # Print night CCE
    # cce_night = abs(carbon_biomass_dict[model.reactions.get_by_id(output_string + "_12").id]) / day_to_night_carbons
    # print("night CCE:", cce_night)

    carbon_flux_dict = {"CO2 exchange": model.reactions.CO2_tx.flux,
                        "Carbon day to night": day_to_night_carbons,
                        "Carbon night to day": night_to_day_carbons,
                        # "CO2 exchange night:day": night_to_day_carbons / day_to_night_carbons,
                        # "total Carbon biomass flux": sum(carbon_biomass_dict.values()),
                        # "total CCE": CCE,
                        # "nCO2RC": nCO2RC,
                        # "night CCE": cce_night
                        }

    # create a pandas dataframe from the dictionary
    carbon_flux_df = pd.DataFrame.from_dict(carbon_flux_dict, orient='index', columns=['value'])

    return carbon_flux_df
    


### calculating how much of nocturnally produced CO2 is refixed and how much is released ###

def calculate_noct_CO2_refixation_coefficient(model, verbose=False, threshold = 0.00001, show_noflux = False):

    CO2_c2_produced = 0
    CO2_c2_consumed = 0
    
    organelles = {"Cytosol":"_c","Plastid": "_p","Mitochondrion": "_m", "Peroxisome":"_x"}

    CO2_rxns = []

    for rxn in model.reactions:
        #find all nocturnal reactions that involve CO2, except the import from external into cytosol
        for met in rxn.metabolites:
            current_met = met.id
            if ("CARBON_DIOXIDE_" in met.id or "HCO3_" in met.id) and met.id[-3:]=="_12" and "CO2_ec" not in rxn.id and rxn.id not in CO2_rxns:
                CO2_rxns.append(rxn.id)
            # if re.search('CARBON_DIOXIDE_._12', met.id) and (rxn.id != "CO2_tx_12" and rxn.id != "CO2_ec_12"):
            if re.search('CARBON_DIOXIDE_._12', met.id) and (rxn.id != "CO2_tx_12" and rxn.id != "CO2_ec_12" and rxn.id != "CO2_mc_12" and rxn.id != "CO2_pc_12" and rxn.id != "CO2_xc_12"):
                CO2_flux = rxn.metabolites.get(model.metabolites.get_by_id(current_met)) * rxn.flux
                if CO2_flux > 0:
                    CO2_c2_produced += CO2_flux
                    # print("produced:", rxn.id, CO2_flux)

                elif CO2_flux < 0:
                    CO2_c2_consumed += CO2_flux
                    # print("consumed:", rxn.id, CO2_flux)

    for item in organelles:
        if verbose:
            print("-----", item, "-----")
        for rxn in CO2_rxns:
            for met in model.reactions.get_by_id(rxn).metabolites:
                if "CARBON_DIOXIDE"+organelles[item] in met.id: 
                    CO2_flux = model.reactions.get_by_id(rxn).flux * model.reactions.get_by_id(rxn).get_coefficient(met)
                    if CO2_flux > 0 and abs(CO2_flux) > threshold:
                        if verbose:
                            print("\033[96m"+"producing CO2:\t"+ rxn+ "\t"+str(round(model.reactions.get_by_id(rxn).flux, 6))+"\033[00m")
                    elif CO2_flux < 0 and abs(CO2_flux) > threshold:
                        if verbose:
                            print("\033[91m"+"consuming CO2:\t"+ rxn+ "\t"+str(round(model.reactions.get_by_id(rxn).flux, 6))+"\033[00m")
                    else:
                        if verbose and show_noflux:
                            print("\033[97m"+"not consuming/producing CO2:\t"+ rxn + "\t"+str(model.reactions.get_by_id(rxn).flux) + "\033[00m")
                elif "HCO3"+organelles[item] in met.id: 
                    CO2_flux = model.reactions.get_by_id(rxn).flux * model.reactions.get_by_id(rxn).get_coefficient(met)
                    if CO2_flux > 0 and abs(CO2_flux) > threshold:
                        if verbose:
                            print("\033[96m"+"producing HCO3-:\t"+ rxn+ "\t"+str(round(model.reactions.get_by_id(rxn).flux, 6))+"\033[00m")
                    elif CO2_flux < 0 and abs(CO2_flux) > threshold:
                        if verbose:
                            print("\033[91m"+"consuming HCO3-:\t"+ rxn+ "\t"+str(round(model.reactions.get_by_id(rxn).flux, 6))+"\033[00m")
                    else:
                        if verbose and show_noflux:
                            print("\033[97m"+"not consuming/producing HCO3-:\t"+ rxn + "\t"+str(model.reactions.get_by_id(rxn).flux) + "\033[00m")
                            
                

    CO2_exchange = model.reactions.CO2_ec_12.flux


    if verbose:
        print("\nCO2 produced:", round(CO2_c2_produced, 5))
        print("CO2 consumed:", round(CO2_c2_consumed,5))
        print("CO2 exchange:", round(CO2_exchange, 5))

    try:
        CO2_metric = abs(CO2_c2_consumed)/abs(CO2_c2_produced)
    except:
        CO2_metric = 0
    if verbose:
        print("\nnocturnal CO2 refixation coeffient:", round(CO2_metric,5))
    return CO2_metric
    

def calculate_CO2_refixation_coefficient(model, verbose=False, threshold = 0.00001, show_noflux = False, time_tag = "_12"):

    """
    function to calculate how much of metabolically produced CO2 is refixed and how much is released
    
    input:
    - model: an optimised model
    - threshold: minimum flux to consider reaction flux high enough to print as consuming/producing
    - show_noflux: if True, prints even names of CO2-involving reactions that don't carry any flux or flux above threshold
    - time_tag: the string added to the end of reactions to represent the time phase (e.g. "_night"
    
    output:
    - CO2_metric: CO2 refixation metric (sum of all CO2 consumed by metabolic reactions in specified phase divided by 
                  sum of all CO2 produced by metabolic reactions in specified phase)
    """
    
    CO2_produced = 0
    CO2_consumed = 0
    
    organelles = {"Cytosol":"_c","Plastid": "_p","Mitochondrion": "_m", "Peroxisome":"_x",  "Endoplasmic reticulum":"_r", "Vacuole":"_v", 
                  "Mitochondrion innermembrane interacting with cristal space":"_mi", "Mitochondrion innermembrane interacting with inter membrane space":"_mc",
                 "Mitochondrial intermembrane space":"_i", "Thylakoid":"l"}

    CO2_rxns = []

    #find all nocturnal reactions that involve CO2, except the import from external into cytosol
    for rxn in model.reactions:
        for met in rxn.metabolites:
            current_met = met.id
            
            # if a reaction involves CO2, add it to the list of CO2 reactions (excluding reaction that releases cytosolic CO2 to external environment)
            if "CARBON_DIOXIDE_" in met.id and met.id.endswith(time_tag) and "CO2_ec" not in rxn.id and rxn.id not in CO2_rxns:
                CO2_rxns.append(rxn.id)

            # if a reaction involves CO2, and is not a CO2 transport reaction between compartments, calculate how much CO2 is consumed or produced by it, and add to corresponding running total
            if re.search('CARBON_DIOXIDE_.'+time_tag, met.id) and (rxn.id != "CO2_tx"+time_tag and rxn.id != "CO2_ec"+time_tag 
                                                                   and rxn.id != "CO2_mc"+time_tag and rxn.id != "CO2_pc"+time_tag 
                                                                   and rxn.id != "CO2_xc"+time_tag and rxn.id != "CARBON_DIOXIDE_rc"+time_tag):
                
                CO2_flux = rxn.metabolites.get(model.metabolites.get_by_id(current_met)) * rxn.flux
                
                if CO2_flux > 0:
                    CO2_produced += CO2_flux
                    # print("produced:", rxn.id, CO2_flux)

                elif CO2_flux < 0:
                    CO2_consumed += CO2_flux
                    # print("consumed:", rxn.id, CO2_flux)
    if verbose:
        for item in organelles:
            print("-----", item, "-----")
            for rxn in CO2_rxns:
                for met in model.reactions.get_by_id(rxn).metabolites:
                    if "CARBON_DIOXIDE"+organelles[item] in met.id: 
                        CO2_flux = model.reactions.get_by_id(rxn).flux * model.reactions.get_by_id(rxn).get_coefficient(met)
                        if CO2_flux > 0 and abs(CO2_flux) > threshold:
                            print("\033[96m"+"producing CO2:\t"+ rxn+ "\t"+str(round(model.reactions.get_by_id(rxn).flux, 6))+"\033[00m")
                        elif CO2_flux < 0 and abs(CO2_flux) > threshold:
                            print("\033[91m"+"consuming CO2:\t"+ rxn+ "\t"+str(round(model.reactions.get_by_id(rxn).flux, 6))+"\033[00m")
                        else:
                            if show_noflux:
                                print("\033[97m"+"not consuming/producing CO2:\t"+ rxn + "\t"+str(model.reactions.get_by_id(rxn).flux) + "\033[00m")
                

                
    CO2_exchange = model.reactions.get_by_id("CO2_ec"+time_tag).flux

    if verbose:
        print("\nCO2 produced:", round(CO2_produced, 5))
        print("CO2 consumed:", round(CO2_consumed,5))
        print("CO2 exchange:", round(CO2_exchange, 5))

        
    #calculate how much CO2
    try:
        CO2_metric = abs(CO2_consumed)/abs(CO2_produced)
    except:
        CO2_metric = 0 #return 0 when no CO2 is produced
        
    if verbose:
        print("\nCO2 refixation coeffient for phase \"" +time_tag+"\": "+ str(CO2_metric)+"  (= ~"+str(round(CO2_metric*100,4))+"%)")
    
    return CO2_metric
    




#########################################

def parameter_scan_CO2(model, objective, process, start, stop, stepsize, pFBA, CO2_refixation_allowed=True, verbose=False, iterations=10, 
                       double_ub = False, m=0.01025582745716949, b =-0.01170011252656307, scan_ub_and_lb = False):
    solution_dict = {}
    
    # make list of all reactions to save in the dataframe
    rxn_list = []
    for rxn in model.reactions:
        rxn_list.append(rxn.id)
    # add other metrics to save in dataframe
    metrics_list = ["noct. CO2 coeff.", 
                    "CCE", 
                    "sum of fluxes", 
                    "photon use", 
                    "growth rate ub", 
                    "photon use ub", 
                    "daytime CO2 uptake ub", 
                    "scan process ub", 
                    "total aa day", 
                    "total aa night",
                    "nighttime CCE"]
    for metric in metrics_list:
        rxn_list.append(metric)
        
    
    # create dataframe (rows: reactions and other metrics; columns: flux/value for each model iteration with specified scan value)
    model_rxns = {"reactions": rxn_list}
    df = pd.DataFrame(model_rxns)
    
    # start scan by changing upper bound of specific reaction to current scan value
    print("Scanning process:", process,"\n")
    
    for count, scan_value in enumerate(np.arange(start,stop,stepsize)):
        temp_model = model.copy()
        
        print("------- Scan iteration:", count+1,"    ", "Scan value:", scan_value, "("+process+") --------")
        
        
        # try to optimise the model
        fluxes_list =[]
        try:
            if CO2_refixation_allowed:
                if double_ub:
                    temp_model.reactions.get_by_id(process).bounds = (0, scan_value)
                    temp_model.reactions.get_by_id(objective).bounds = (0, scan_value*m +b)

                    solution_temp, opt_model = optimise_model(temp_model, objective, pFBA)
                else:
                    if scan_ub_and_lb:
                        temp_model.reactions.get_by_id(process).bounds = (scan_value, scan_value)
                    else:
                        temp_model.reactions.get_by_id(process).bounds = (0, scan_value)

                    solution_temp, opt_model = optimise_model(temp_model, objective, pFBA)
            else:
                solution_temp, opt_model, _ = estimateOutputFromNetCO2_Corinna_dev(temp_model,netCO2uptake=scan_value,Output_ID=objective,
                                      Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_00",
                                      CO2in_ID="CO2_tx_00",verbose=verbose, iterations=iterations, pFBA=pFBA)
        except:
            # if infeasible, save every flux and metric as NaN to add to dataframe
            for rxn in temp_model.reactions:
                fluxes_list.append(np.nan)
            for metric in metrics_list:
                fluxes_list.append(np.nan)
            print("\033[91m"+"Model iteration is infeasible!\033[00m")
            
            
        else:
        # runs when model optimisation produced feasible result
        # saves fluxes of each reactions, and value of metric to add to dataframe
            print("Model iteration feasible. Objective flux:", round(opt_model.reactions.get_by_id(objective).flux,7))
            for rxn in opt_model.reactions:
                flux = rxn.flux
                if "linker_12_to_00" in rxn.id:
                    flux = (-1)*flux #to make night to day accumulation negative for plotting later

                fluxes_list.append(flux)
            
            #noct. CO2 coeff.
            CO2_coeff = calculate_noct_CO2_refixation_coefficient(opt_model, verbose=verbose)
            fluxes_list.append(CO2_coeff)
            
            #CCE
            CCE = calculate_CCE(opt_model)
            fluxes_list.append(CCE)
            
            #calculating sum of all fluxes
            sum_of_fluxes = 0
            for item in solution_temp:
                sum_of_fluxes += abs(item)
            fluxes_list.append(sum_of_fluxes)
            
            #calculating photon use
            fluxes_list.append(opt_model.reactions.Photon_tx_00.flux/opt_model.reactions.Photon_tx_00.upper_bound)
            
            #objective, photon, CO2 upper bound
            fluxes_list.append(opt_model.reactions.get_by_id(objective).upper_bound)
            fluxes_list.append(opt_model.reactions.Photon_tx_00.upper_bound)
            fluxes_list.append(opt_model.reactions.CO2_tx_00.upper_bound)
            
            #scan reaction upper bound
            fluxes_list.append(opt_model.reactions.get_by_id(process).upper_bound)
            
            #total amino acid accumulation
            amino_acids = ["GLN_v","VAL_v","PHE_v","TYR_v","L_ALPHA_ALANINE_v","SER_v",
                       "TRP_v","GLT_v","ILE_v","L_ASPARTATE_v","CYS_v","MET_v","bHIS_v",
                       "ARG_v","THR_v","LEU_v","PRO_v","LYS_v","ASN_v","GLY_v"]
            aa_sum_day = 0
            aa_sum_night = 0
            for item in amino_acids:
                for rxn in opt_model.reactions:
                    if item in rxn.id and "linker_00_to_12" in rxn.id:
                        aa_sum_day += rxn.flux
                    elif item in rxn.id and "linker_12_to_00" in rxn.id:
                        aa_sum_night += rxn.flux
            fluxes_list.append(aa_sum_day)
            fluxes_list.append((-1)*aa_sum_night)
                        

            solution_dict[count]= solution_temp

            #nighttime CCE
            day_to_night_carbons  = sum(rxn.flux * rxn.reactants[0].elements.get("C", 0) for rxn in opt_model.reactions if "linker_00_to_12" in rxn.id)
            night_CCE = 1 + (opt_model.reactions.CO2_tx_12.flux / day_to_night_carbons)
            fluxes_list.append(night_CCE)
        
        #name model iteration after current scan value
        column_name = str(round(scan_value, 3))
        df.insert(count+1, column_name, fluxes_list)

    
    df.set_index("reactions", inplace= True)
    
    print("\n")
    return df, solution_dict

def parameter_scan_CO2_CAM(model, objective, process, start, stop, stepsize, pFBA, CO2_refixation_allowed=True, 
                           verbose=False, iterations=10, double_ub = False, CAM_constraint=False):
    solution_dict = {}
    
    # make list of all reactions to save in the dataframe
    rxn_list = []
    for rxn in model.reactions:
        rxn_list.append(rxn.id)
    # add other metrics to save in dataframe
    metrics_list = ["noct. CO2 coeff.", "CCE", "sum of fluxes", "photon use", "growth rate ub", "photon use ub", "daytime CO2 uptake ub", "scan process ub", "total aa day", "total aa night"]
    for metric in metrics_list:
        rxn_list.append(metric)
        
    
    # create dataframe (rows: reactions and other metrics; columns: flux/value for each model iteration with specified scan value)
    model_rxns = {"reactions": rxn_list}
    df = pd.DataFrame(model_rxns)
    
    # start scan by changing upper bound of specific reaction to current scan value
    print("Scanning process:", process,"\n")
    
    for count, scan_value in enumerate(np.arange(start,stop,stepsize)):
        temp_model = model.copy()
        
        print("------- Scan iteration:", count+1,"    ", "Scan value:", scan_value, "("+process+") --------")
        
        
        # try to optimise the model
        fluxes_list =[]
        try:
            if CO2_refixation_allowed:
                if double_ub:
                    temp_model.reactions.get_by_id(process).bounds = (0, scan_value)
                    temp_model.reactions.get_by_id(objective).bounds = (0, scan_value*0.0096707943407801 -0.0125515499006004)

                    solution_temp, opt_model = optimise_model(temp_model, objective, pFBA)
                elif CAM_constraint:
                    temp_model.reactions.get_by_id(process).bounds = (0, scan_value)
                    
                    # temp_model.reactions.get_by_id("PEPCARBOX_RXN_c_12").bounds = (0, scan_value*2.094050238268857 -0.004094236753110575)
                    
                    try:
                        model.remove_cons_vars(model.constraints.PEPC_obj_const)
                        print("PEPC_obj_const constraint removed")
                    except:
                        print("")        

                    PEPC_obj_const = model.problem.Constraint(
                            (2.094050238268857*(model.reactions.get_by_id(objective).flux_expression)-0.004094236753110575) -
                            (model.reactions.PEPCARBOX_RXN_c_12.flux_expression) , 
                            lb = 0,
                            ub= 1000,
                            name= "PEPC_obj_const")
                    model.add_cons_vars(PEPC_obj_const)
                    
                    temp_model.reactions.get_by_id("PEPCARBOX_RXN_p_12").bounds = (0, 0)

                    
                    solution_temp, opt_model = optimise_model(temp_model, objective, pFBA)
                else:
                
                    temp_model.reactions.get_by_id(process).bounds = (0, scan_value)

                    solution_temp, opt_model = optimise_model(temp_model, objective, pFBA)
            else:
                solution_temp, opt_model, _ = estimateOutputFromNetCO2_Corinna_dev(temp_model,netCO2uptake=scan_value,Output_ID=objective,
                                      Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_00",
                                      CO2in_ID="CO2_tx_00",verbose=verbose, iterations=iterations, pFBA=pFBA)
        except:
            # if infeasible, save every flux and metric as NaN to add to dataframe
            for rxn in temp_model.reactions:
                fluxes_list.append(np.nan)
            for metric in metrics_list:
                fluxes_list.append(np.nan)
            print("\033[91m"+"Model iteration is infeasible!\033[00m")
            
            
        else:
        # runs when model optimisation produced feasible result
        # saves fluxes of each reactions, and value of metric to add to dataframe
            print("Model iteration feasible. Objective flux:", round(opt_model.reactions.get_by_id(objective).flux,7))
            for rxn in opt_model.reactions:
                flux = rxn.flux
                if "linker_12_to_00" in rxn.id:
                    flux = (-1)*flux #to make night to day accumulation negative for plotting later

                fluxes_list.append(flux)
            
            #noct. CO2 coeff.
            CO2_coeff = calculate_noct_CO2_refixation_coefficient(opt_model, verbose=verbose)
            fluxes_list.append(CO2_coeff)
            
            #CCE
            CCE = calculate_CCE(opt_model)
            fluxes_list.append(CCE)
            
            #calculating sum of all fluxes
            sum_of_fluxes = 0
            for item in solution_temp:
                sum_of_fluxes += abs(item)
            fluxes_list.append(sum_of_fluxes)
            
            #calculating photon use
            fluxes_list.append(opt_model.reactions.Photon_tx_00.flux/opt_model.reactions.Photon_tx_00.upper_bound)
            
            #objective, photon, CO2 upper bound
            fluxes_list.append(opt_model.reactions.get_by_id(objective).upper_bound)
            fluxes_list.append(opt_model.reactions.Photon_tx_00.upper_bound)
            fluxes_list.append(opt_model.reactions.CO2_tx_00.upper_bound)
            
            #scan reaction upper bound
            fluxes_list.append(opt_model.reactions.get_by_id(process).upper_bound)
            
            #total amino acid accumulation
            amino_acids = ["GLN_v","VAL_v","PHE_v","TYR_v","L_ALPHA_ALANINE_v","SER_v",
                       "TRP_v","GLT_v","ILE_v","L_ASPARTATE_v","CYS_v","MET_v","bHIS_v",
                       "ARG_v","THR_v","LEU_v","PRO_v","LYS_v","ASN_v","GLY_v"]
            aa_sum_day = 0
            aa_sum_night = 0
            for item in amino_acids:
                for rxn in opt_model.reactions:
                    if item in rxn.id and "linker_00_to_12" in rxn.id:
                        aa_sum_day += rxn.flux
                    elif item in rxn.id and "linker_12_to_00" in rxn.id:
                        aa_sum_night += rxn.flux
            fluxes_list.append(aa_sum_day)
            fluxes_list.append((-1)*aa_sum_night)
                        

            solution_dict[count]= solution_temp
        
        #name model iteration after current scan value
        column_name = str(round(scan_value, 3))
        df.insert(count+1, column_name, fluxes_list)

    
    df.set_index("reactions", inplace= True)
    
    print("\n")
    return df, solution_dict

def set_PEPC_C3_constraints(model, slope=2.1844080596948623, objective = "Marchantia_biomass"):
    """
    Sets the ratio between the objective reaction flux and the flux through PEP_carboxylase_c at night
    """
    try:
        model.remove_cons_vars(model.constraints.PEPC_obj_const)
        print("previous PEPC_obj_const constraint removed")
    except:
        print("")        

    PEPC_obj_const = model.problem.Constraint(
            (slope*(model.reactions.get_by_id(objective).flux_expression)) -
            (model.reactions.PEPCARBOX_RXN_c_12.flux_expression) , 
            lb = 0,
            ub= 1000, #ub allows PEPC_c_night flux to be smaller than in C3 type scenario, but not higher
            name= "PEPC_obj_const")
    model.add_cons_vars(PEPC_obj_const)
    print("Setting new PEPC_obj_const constraint, with slope", str(slope),".")
    
    model.reactions.get_by_id("PEPCARBOX_RXN_p_12").bounds = (0, 0)
                          
    return slope, ""






def optimise_model(model, objective_reaction, pFBA=True, direction="max"):
    """
    Optimize the given model for the specified objective reaction using pFBA or FBA.

    Args:
        model (cobra.Model): The metabolic model to optimize.
        objective_reaction (str): The ID of the reaction to use as the objective.
        pFBA (bool, optional): Whether to use pFBA (True) or FBA (False). Defaults to True.

    Returns:
        tuple: A tuple containing the solution and the updated model.

    Raises:
        AttributeError: If the objective_reaction is not present in the model.
    """
    print("Objective:", objective_reaction)
    
    try:
        reaction = model.reactions.get_by_id(objective_reaction)
    except AttributeError:
        raise AttributeError(f"Reaction {objective_reaction} not found in model.")
    
    model.objective = {reaction: 1}
    model.objective_direction = direction
    if pFBA:
        solution = cobra.flux_analysis.pfba(model)
    else:
        solution = model.optimize()
    return solution, model


def extract_fluxes_from_linker_reactions(dataframe, threshold=0.00001):
    """
    Extracts fluxes from linker reactions in a pandas DataFrame and returns them as dictionaries.
    """
    y_values = {}
    y_values_pos = {}
    y_values_neg = {}

    for index, item in dataframe.iterrows():
        if "linker" in index:
            # sum_flux = item.sum()
            # if abs(sum_flux) > threshold:
            max_flux = max(abs(item.max()),abs(item.min()))
            if abs(max_flux) > threshold:
                y_values[index] = item.round(5).tolist()
                y_values_pos[index] = item.clip(lower=0).round(5).tolist()
                y_values_neg[index] = item.clip(upper=0).round(5).tolist()

    return y_values, y_values_pos, y_values_neg

def plot_accum(model, objective_reaction, dataframe, threshold, growth_rate, CO2_day_rate, lim2_x=0, lim2_y=1.4, lim3_x=-1.5, lim3_y=7.2, lim5_x=-0.05, lim5_y=0.05, time_phase="_12", CO2_organelle = "_c"):
    results = dataframe.columns.tolist()
    threshold = 0.0001

    ### Extracting fluxes from linker reactions ##################################################################################
    y_values, y_values_pos, y_values_neg = extract_fluxes_from_linker_reactions(dataframe, threshold)

    ### Extracting reactions that involve nocturnal, cytosolic CO2 ##################################################################################

    CO2_12_fluxes = {}
    total_CO2_produced = [0] * len(dataframe.columns)
    total_CO2_consumed = [0] * len(dataframe.columns)

    CO2_metabolite = "CARBON_DIOXIDE"

    # iterate over each row in the dataframe
    for index, item in dataframe.iterrows():
        CO2_x_list = []
        try:
            # iterate over reactions in the model to find any that contain the CO2 metabolite
            for met in model.reactions.get_by_id(index).metabolites:
                if met.id[:-5] == CO2_metabolite and met.id[-5:-3] in CO2_organelle and met.id[-3:] == time_phase and "_pc_" not in index and "_xc_" not in index and "_mc_" not in index and "_ec_" not in index:
                    CO2_met = met.id
                    sum = 0
                    
                    # for each CO2 metabolite containing reaction, get the flux for each model iteration and sum them up separately for production and consumption
                    for count, x in enumerate(item):
                        if not np.isnan(x):
                            y = x * model.reactions.get_by_id(index).get_coefficient(CO2_met)
                            sum += y
                            CO2_x_list.append(y)   
                            if y > 0:
                                total_CO2_produced[count] += y
                            elif y < 0:
                                total_CO2_consumed[count] += y
                        else:
                            CO2_x_list.append(np.nan)   
                    CO2_12_fluxes[index] = CO2_x_list

            CO2_12_fluxes["total CO2 produced"] = total_CO2_produced
            CO2_12_fluxes["total CO2 consumed"] = total_CO2_consumed
        except:
            #exception needed in cases where metric was added to dataframe (eg noct CO2 refixation coefficient), as they don't have .metabolites attribute
            continue
            



    ######################## Plotting ###########################################################################################################################           

    fig, (ax2, ax3, ax5, ax4) = plt.subplots(4, sharex = True)

    ax2.set_title("Objective")
    ax2.plot(results, dataframe.loc[objective_reaction], label="growth rate", linewidth=2.2)

    ax2.plot(results, dataframe.loc["sum of fluxes"]/1000, label="sum of fluxes/1000")
    ax2.axhspan(growth_rate[0], growth_rate[1], color='gray', alpha=0.3, lw=0)
    ax2.set_ylim((0,1.1))
    ax2.plot(results, dataframe.loc["photon use"], label="photon use")
    ax2.plot(results, dataframe.loc["growth rate ub"], "--", label="growth rate ub", linewidth=1.2, color="gray")

    ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))


    ax3.set_title("CO2 Exchange")
    ax3.plot(results, dataframe.loc["CO2_tx_00"], label="CO2 day", linewidth=2)
    ax3.plot(results, dataframe.loc["daytime CO2 uptake ub"], "--", label="CO2 day ub", linewidth=1.2, color="gray")

    ax3.plot(results, dataframe.loc["CO2_tx_12"], label="CO2 night")
    ax3.plot(results, dataframe.loc["noct. CO2 coeff."], label="CO2 refix. coeff.", linewidth=2.2)
    ax3.plot(results, dataframe.loc["CCE"],"-.", label="CCE")
    ax3.plot(results, dataframe.loc["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_00"], label="Rubisco")
    ax3.plot(results, dataframe.loc["ISOCITDEH_RXN_c_12"], label="Isocit deh_c night")
    ax3.plot(results, dataframe.loc["RXN0_5224_c_12"], label="carb anhyd_c night")
    ax3.plot(results, dataframe.loc["PEPCARBOX_RXN_c_12"],"--", label="PEPC_c night")
    ax3.plot(results, dataframe.loc["PEPCARBOX_RXN_c_00"],"--", label="PEPC_c day")

    # ax3.plot(results, dfx.loc["PEPCARBOXYKIN_RXN_c_12"], label="PEPCkinase_c night")

    ax3.axhspan(CO2_day_rate[0], CO2_day_rate[1], color='gray', alpha=0.3, lw=0)

    ax3.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    

    #display nocturnal CO2 involving reaction fluxes
    ax4.set_title("nocturnal CO2 production/consumption")
    
    LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
    NUM_COLORS = len(CO2_12_fluxes)
    cm = plt.get_cmap('tab20')
    ax4.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)], linestyle= [LINE_STYLES[i%len(LINE_STYLES)] for i in range(NUM_COLORS)])

    for rxn in CO2_12_fluxes:
        width = 2
        sum_CO2 = 0
        for item in CO2_12_fluxes[rxn]:
            sum_CO2 += item
        if sum_CO2 == 0:
            width = 0.05
        if sum_CO2 != 0:
            if "total" in rxn:
                ax4.plot(results, CO2_12_fluxes[rxn],"-.", label=rxn, linewidth=1)
            else:
                ax4.plot(results, CO2_12_fluxes[rxn], label=rxn, linewidth=width)
    ax4.legend(loc="center left", bbox_to_anchor=(1, 0.5))


    # display accumulation related fluxes
    ax5.set_title("Accumulation")

    LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
    NUM_STYLES = len(LINE_STYLES)
    NUM_COLORS = len(y_values)
    ax5.set_title("Accumulation")
    cm = plt.get_cmap('gist_rainbow')
    ax5.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)], linestyle= [LINE_STYLES[i%len(LINE_STYLES)] for i in range(NUM_COLORS)])
    ax5.plot(results, dataframe.loc["total aa day"],"--", label="total aa day", linewidth=1.2, color="gray")
    ax5.plot(results, dataframe.loc["total aa night"],"--", label="total aa night", linewidth=1.2, color="gray")

    


    for values in y_values:
        ax5.plot(results, y_values[values], label=values[:-16], linewidth=2)
    ax5.legend(loc="center left", bbox_to_anchor=(1, 0.5))


    ax3.set_xticks(range(len(dataframe.columns)))
    ax3.set_xticklabels(dataframe.columns, rotation=90)
    ax2.set_xticks(range(len(dataframe.columns)))
    ax2.set_xticklabels(dataframe.columns, rotation=90, alpha=0)

    plt.xticks(range(len(dataframe.columns)), rotation=90, label=dataframe.columns)

    fig.set_figwidth(10)
    fig.set_figheight(24)
    
    ax2.set_ylim((lim2_x, lim2_y))
    ax3.set_ylim((lim3_x,lim3_y))
    ax5.set_ylim((lim5_x,lim5_y))







import matplotlib.pyplot as plt
import numpy as np

def plotBoundary(Model, solutions,thresh=1e-5,displayModelID=0,solutionNames=[]):
    bdryRxns = sorted([x for x in Model.boundary if any([abs(sol[x.id])>thresh for sol in solutions])],key=lambda x:abs(solutions[0][x.id]),reverse=True)
    width = 0.9/len(solutions)
    plt.ylabel('flux')
    for it,sol in enumerate(solutions):
        plt.bar(np.arange(len(bdryRxns))+(it)*width,[sol[x.id] for x in bdryRxns],width=width)
    if displayModelID:
        plt.xticks(range(len(bdryRxns)),[x.id for x in bdryRxns],rotation=90)
    else:
        plt.xticks(range(len(bdryRxns)),[x.name for x in bdryRxns],rotation=90)
    if solutionNames and len(solutionNames)==len(solutions):
        plt.legend(solutionNames)
    else:
        if solutionNames and len(solutionNames)!=len(solutions):
            print('Length of solutionNames does not match length of solutions.')
        plt.legend(['Solution '+str(x) for x in range(len(solutions))])

def plotATPBudget(Model,solutions,thresh=1e-5,displayModelID=0,solutionNames=[]):
    ATP_rxns=sorted([x for x in Model.metabolites.atp_c.reactions if any([abs(sol[x.id])>thresh for sol in solutions])],key=lambda x: abs(solutions[0][x.id]*x.get_coefficient('ATP_c')),reverse=True)
    cm = plt.get_cmap('tab20')
    posSum=np.zeros(len(solutions))
    negSum=np.zeros(len(solutions))
    width = 0.7

    thisplot = plt.figure(figsize=(3,4))
    plts=[]
    for it,rxn in enumerate(ATP_rxns):
        rxnFluxes=[sol[rxn.id]*rxn.get_coefficient('ATP_c') for sol in solutions]
        plts+=[plt.bar(np.arange(len(solutions)), rxnFluxes, bottom=[negSum[x] if rxnFluxes[x]<0 else posSum[x] for x in range(len(rxnFluxes))],color=cm(it/len(ATP_rxns)),width=width)]
        negSum=[negSum[x]+rxnFluxes[x] if rxnFluxes[x]<0 else negSum[x] for x in range(len(negSum))]
        posSum=[posSum[x]+rxnFluxes[x] if rxnFluxes[x]>0 else posSum[x] for x in range(len(posSum))]

    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.xlim(-0.6,len(solutions)-0.4)
    plt.ylabel("ATP produced/consumed (µmol/s)")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    if displayModelID:
        plt.legend(plts,[x.id for x in ATP_rxns],bbox_to_anchor=(1,1))
    else:
        plt.legend(plts,[x.name for x in ATP_rxns],bbox_to_anchor=(1,1))
    if solutionNames and len(solutionNames)==len(solutions):
        plt.xticks(range(len(solutionNames)),[x for x in solutionNames],rotation=90)
    else:
        if solutionNames and len(solutionNames)!=len(solutions):
            print('Length of solutionNames does not match length of solutions.')
    plt.show()    



#### wrong metabolite ids for NAD(P)/H
def plotNADBudget(Model,solutions,thresh=1e-5,displayModelID=0,solutionNames=[]):
    NADH_rxns = list(Model.metabolites.nadh_c.reactions)
    NADPH_rxns = list(Model.metabolites.nadph_c.reactions)
    nadh = Model.metabolites.nadh_c
    nadph = Model.metabolites.nadph_c
    NAD_rxns=sorted([x for x in NADH_rxns+NADPH_rxns if any([abs(sol[x.id])>thresh for sol in solutions])],key=lambda x: [abs(solutions[0][x.id]*x.get_coefficient(nadh)) if nadh in x.metabolites else abs(solutions[0][x.id]*x.get_coefficient(nadph))],reverse=True)
    cm = plt.get_cmap('tab20')
    
    posSum=np.zeros(len(solutions))
    negSum=np.zeros(len(solutions))
    width = 0.7

    thisplot = plt.figure(figsize=(3,4))
    plts=[]
    for it,rxn in enumerate(NAD_rxns):
        nadhSol=np.zeros(len(solutions))
        nadphSol=np.zeros(len(solutions))
        if nadh in rxn.metabolites:
            nadhSol=[sol[rxn.id]*rxn.get_coefficient('nadh_c') for sol in solutions]
        if nadph in rxn.metabolites:
            nadphSol=[sol[rxn.id]*rxn.get_coefficient('nadph_c') for sol in solutions]
        rxnFluxes=[nadhSol[x]+nadphSol[x] for x in range(len(solutions))]
        plts+=[plt.bar(np.arange(len(solutions)), rxnFluxes, bottom=[negSum[x] if rxnFluxes[x]<0 else posSum[x] for x in range(len(rxnFluxes))],color=cm(it/len(NAD_rxns)),width=width)]
        negSum=[negSum[x]+rxnFluxes[x] if rxnFluxes[x]<0 else negSum[x] for x in range(len(negSum))]
        posSum=[posSum[x]+rxnFluxes[x] if rxnFluxes[x]>0 else posSum[x] for x in range(len(posSum))]

    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.xlim(-0.6,len(solutions)-0.4)
    plt.ylabel("NAD(P)H produced/consumed (µmol/s)")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    if displayModelID:
        plt.legend(plts,[x.id for x in NAD_rxns],bbox_to_anchor=(1,1))
    else:
        plt.legend(plts,[x.name for x in NAD_rxns],bbox_to_anchor=(1,1))
    if solutionNames and len(solutionNames)==len(solutions):
        plt.xticks(range(len(solutionNames)),[x for x in solutionNames],rotation=90)
    else:
        if solutionNames and len(solutionNames)!=len(solutions):
            print('Length of solutionNames does not match length of solutions.')
    plt.show()


def add_pathway_to_compartments(plant_model, pathway_model, time_phases, compartments, transporters, linkers,
                                charge_states, charge_states_dict):
    
    pathway_model_comp = cobra.Model()
    
    
    #gets list of all metabolites by ids, in any compartment, that are in the plant model
    plant_mets = []
    for plant_met in plant_model.metabolites:
        if plant_met.id[:-5] not in plant_mets:
            plant_mets.append(plant_met.id[:-5])
            
    #gets list of all ractions by ids, in any compartment, that are in the plant model
    plant_rxns = []
    for plant_rxn in plant_model.reactions:
        if plant_rxn.id[:-5] not in plant_rxns:
            plant_rxns.append(plant_rxn.id[:-5])
    
    
    
    # create metabolites and reactions for specified time phases and compartments
    for time in time_phases:
        for comp in compartments:
            for met in pathway_model.metabolites:
                temp_met = Metabolite(met.id+comp+time, name=met.id+comp+time, compartment=comp[1], charge = met.charge, formula = met.formula)
                pathway_model_comp.add_metabolites(temp_met)
            for rxn in pathway_model.reactions:
                temp_rxn = Reaction(rxn.id+comp+time)
                temp_rxn.upper_bound = rxn.upper_bound
                temp_rxn.lower_bound = rxn.lower_bound
                temp_rxn.subsystem = rxn.subsystem
                temp_rxn.annotation = rxn.annotation
                temp_rxn.notes = rxn.notes
                
                for met in rxn.metabolites:
                    temp_rxn.add_metabolites({pathway_model_comp.metabolites.get_by_id(met.id+comp+time): rxn.get_coefficient(met)})
                if temp_rxn.id not in plant_rxns:
                    pathway_model_comp.add_reactions([temp_rxn])
    
    
    
    # add transporters between cytosol and specified compartments
    if transporters and len(compartments)>1:

        for met in pathway_model_comp.metabolites:
            #if metabolite is not a plant native metabolite and is not in the cytosol, add a transporter reaction from the cytosol to its compartment
            if met.id[:-5] not in plant_mets and met.id[-5:-3] != "_c":
                
                transp_rxn = Reaction(met.id[:-5]+"_c"+met.id[-4:])
                transp_rxn.upper_bound = 1000.0
                transp_rxn.lower_bound = -1000.0
                transp_rxn.add_metabolites({pathway_model_comp.metabolites.get_by_id(met.id[:-5]+"_c"+met.id[-3:]): -1, 
                                            pathway_model_comp.metabolites.get_by_id(met.id): 1})
                pathway_model_comp.add_reactions([transp_rxn])

            # if met.id[:-5] not in plant_mets:
            #     print(met.id, "already in the plant model, not added again")
            #     print(plant_model.metabolites.get_by_id(met.id).name)           #### not working!!!
    
    # add linker reactions for non-native pathway metabolites
    if linkers:
    
        #add transporter to and from the vacuole
        for met in pathway_model.metabolites:
            if met.id not in plant_mets:
                #create vacuolar metabolite
                temp_met = Metabolite(met.id+"_v_00", name=met.id+"_v_00", compartment="v", charge = met.charge, formula = met.formula)
                pathway_model_comp.add_metabolites(temp_met)
                temp_met = Metabolite(met.id+"_v_12", name=met.id+"_v_12", compartment="v", charge = met.charge, formula = met.formula)
                pathway_model_comp.add_metabolites(temp_met)
                
                #create cytosol metabolite
                temp_met = Metabolite(met.id+"_c_00", name=met.id+"_c_00", compartment="c", charge = met.charge, formula = met.formula)
                pathway_model_comp.add_metabolites(temp_met)
                temp_met = Metabolite(met.id+"_c_12", name=met.id+"_c_12", compartment="c", charge = met.charge, formula = met.formula)
                pathway_model_comp.add_metabolites(temp_met)
                
                vac_transp_rxn = Reaction(met.id+"_cv_00")
                vac_transp_rxn.upper_bound = 1000.0
                vac_transp_rxn.lower_bound = -1000.0
                vac_transp_rxn.add_metabolites({pathway_model_comp.metabolites.get_by_id(met.id+"_c_00"): -1, 
                                            pathway_model_comp.metabolites.get_by_id(met.id+"_v_00"): 1})
                pathway_model_comp.add_reactions([vac_transp_rxn])
                
                vac_transp_rxn = Reaction(met.id+"_cv_12")
                vac_transp_rxn.upper_bound = 1000.0
                vac_transp_rxn.lower_bound = -1000.0
                vac_transp_rxn.add_metabolites({pathway_model_comp.metabolites.get_by_id(met.id+"_c_12"): -1, 
                                            pathway_model_comp.metabolites.get_by_id(met.id+"_v_12"): 1})
                pathway_model_comp.add_reactions([vac_transp_rxn])
            
        
        
        #add linker day to night and night to day
        for met in pathway_model_comp.metabolites:
            #if metabolite is not a plant native metabolite and is not in the cytosol, add a transporter reaction from the cytosol to its compartment
            if met.id[:-5] not in plant_mets and met.id[-5:-3] == "_v":
                
                if met.id[-2:] == "00":
                    linker_rxn_day_to_night = Reaction(met.id[:-3]+"_linker_00_to_12")
                    linker_rxn_day_to_night.upper_bound = 1000.0
                    linker_rxn_day_to_night.lower_bound = 0
                    linker_rxn_day_to_night.add_metabolites({pathway_model_comp.metabolites.get_by_id(met.id):-1, 
                                                pathway_model_comp.metabolites.get_by_id(met.id[:-2]+"12"): 1})
                    pathway_model_comp.add_reactions([linker_rxn_day_to_night])
                elif met.id[-2:] == "12":
                    linker_rxn_night_to_day = Reaction(met.id[:-3]+"_linker_12_to_00")
                    linker_rxn_night_to_day.upper_bound = 1000.0
                    linker_rxn_night_to_day.lower_bound = 0
                    linker_rxn_night_to_day.add_metabolites({pathway_model_comp.metabolites.get_by_id(met.id):-1, 
                                                pathway_model_comp.metabolites.get_by_id(met.id[:-2]+"00"): 1})
                    pathway_model_comp.add_reactions([linker_rxn_night_to_day])
    
    
    # modify metabolites according to defined charge states for different compartments
    if charge_states:

        #add charge state metabolites to the model
        for item in charge_states_dict:
            met = Metabolite("a"+item)
            met.name = "a"+item
            met.compartment = item[-1]
            pathway_model_comp.add_metabolites(met)
            
        for (comp, time) in zip(compartments, time_phases):
            met = Metabolite("PROTON"+comp+time)
            met.name = "PROTON"+comp+time
            met.compartment = comp[-1]
            pathway_model_comp.add_metabolites(met)

        for rxn in pathway_model_comp.reactions:
            charge_state = False
            proton_coeff = 0
            for met in rxn.metabolites:
                # print(met.id)
                if met.id in charge_states_dict:
                    charge_state = True
                    temp_coeff = rxn.get_coefficient(met.id)
                    amet = pathway_model_comp.metabolites.get_by_id("a"+met.id)
                    rxn.subtract_metabolites({met:rxn.get_coefficient(met.id)})

                    rxn.add_metabolites({met:charge_states_dict[met.id]*temp_coeff})
                    rxn.add_metabolites({amet:round((1-charge_states_dict[met.id]),2)*temp_coeff})
                    proton_coeff += (1-charge_states_dict[met.id])*temp_coeff
            if charge_state:
                rxn.add_metabolites({pathway_model_comp.metabolites.get_by_id("PROTON"+met.id[-5:]):-round(proton_coeff,2)})

    
    
    # merge the compartmentalised pathway model with the plant model
    plant_with_pathway_comp_model = plant_model.merge(pathway_model_comp, inplace = False)
    
    
    return plant_with_pathway_comp_model, pathway_model_comp



#####################################################################
def generateCO2budget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="1",
        save_plot_to="temp.png",threshold = 0.01, colourDict={}):
  
  if outfile!="":
    fout = open(outfile,"w")
  CO2dict = dict()
  total = 0
  for p in ("c","p","m","x","e"):
    met=model.metabolites.get_by_id("CARBON_DIOXIDE_"+p+day_or_night_tag)
    for rxn in met.reactions:
      if rxn.id.__contains__("CO2_mc") or rxn.id.__contains__("CO2_pc") or rxn.id.__contains__("CO2_ec") or rxn.id.__contains__("CO2_xc"):
        continue
      sto=rxn.metabolites.get(met)
      if outfile!="":
        fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution[rxn.id]*(sto))+"\t"+met.compartment+"\n")
      CO2dict[rxn.id]=solution[rxn.id]*(sto)
      if solution[rxn.id]*(sto) > 0:
        total = total + (solution[rxn.id]*(sto))
  print("Total:", total)
  if outfile!="":
    fout.close()

  tempDict = dict()
  for rxn in CO2dict.keys():
    tempDict[rxn]=abs(CO2dict[rxn])

  #sort CO2dict by values
  import operator
  sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)

  CO2dict2 = dict()
  CO2dict2["Others-pos"]=0
  CO2dict2["Others-neg"]=0
  baseline = dict()
  pos_base=0
  neg_base=0
  i=0
  for TEMP in sorted_by_value:
    rxn = TEMP[0]
    if CO2dict[rxn]>0:
      if CO2dict[rxn] < total*threshold:
        if percentage:
          CO2dict2["Others-pos"]=CO2dict2["Others-pos"]+float(CO2dict[rxn]*100)/total
        else:
          CO2dict2["Others-pos"]=CO2dict2["Others-pos"]+CO2dict[rxn]
        continue
      base = pos_base
      if percentage:
        CO2dict2[rxn]=float(CO2dict[rxn]*100)/total
        pos_base = pos_base + float(CO2dict[rxn]*100)/total
      else:
        pos_base = pos_base + CO2dict[rxn]
        CO2dict2[rxn]=CO2dict[rxn]
    else:
      if abs(CO2dict[rxn]) < total*threshold:
        if percentage:
          CO2dict2["Others-neg"]=CO2dict2["Others-neg"]+float(CO2dict[rxn]*100)/total
        else:
          CO2dict2["Others-neg"]=CO2dict2["Others-neg"]+CO2dict[rxn]
        continue
      base = neg_base
      if percentage:
        CO2dict2[rxn]=float(CO2dict[rxn]*100)/total
        neg_base = neg_base + float(CO2dict[rxn]*100)/total
      else:
        neg_base = neg_base + CO2dict[rxn]
        CO2dict2[rxn]=CO2dict[rxn]
    i=i+1
    baseline[rxn]=base
  baseline["Others-pos"]=pos_base
  baseline["Others-neg"]=neg_base

  if show_plot:
    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.figure(figsize=(3,4))
    for rxn in CO2dict2.keys():
      if colourDict.keys().__contains__(rxn):
        plt.bar(1,CO2dict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn])
      else:
        plt.bar(1,CO2dict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
    plt.xlim(0.8,1.2)
    if percentage:
      plt.ylabel("CO2 produced/consumed (%)")
    else:
      plt.ylabel("CO2 produced/consumed (in moles)")
    handles, labels = plt.gca().get_legend_handles_labels()
    labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
    handles2=[handles[labels.index(i)] for i in labels2]
    lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')

  # create a pandas dataframe from the dictionary
  CO2_df = pd.DataFrame.from_dict(CO2dict2, orient='index', columns=['value'])

  return CO2dict2, CO2_df

#####################################################################
def generateMetaboliteBudget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="1",
        save_plot_to="temp.png",threshold = 0.02, colourDict={}, metabolite=" ", compartments=("c","p","m","x","e")):
  
  if outfile!="":
    fout = open(outfile,"w")
  metdict = dict()
  
  transport_rxns = ("_pc_", "_mc_", "_ec_", "_xc_")

  total = 0

  for p in compartments:
    try:
        met=model.metabolites.get_by_id(metabolite+p+day_or_night_tag)
    except:
        # print(met +"not present in"+ p)
        pass
    else:
        for rxn in met.reactions:

            sto=rxn.metabolites.get(met)

            if outfile!="":
                fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution[rxn.id]*(sto))+"\t"+met.compartment+"\n")
            
            #if reactions is not already in the dictionary, add it
            #if reaction is in dictionary, add the flux times the metabolite stoichiometry to the existing value
                # this is to account for metabolite transport reactions
            if rxn.id not in metdict.keys():
              metdict[rxn.id]=solution[rxn.id]*(sto)
            else:
              metdict[rxn.id]=metdict[rxn.id]+solution[rxn.id]*(sto)
            
            #add the flux times the metabolite stoichiometry to the total
            if solution[rxn.id]*(sto) > 0 and not any(c in rxn.id for c in transport_rxns):
                total = total + (solution[rxn.id]*(sto))
                
  print("Total:", total)

  if outfile!="":
    fout.close()

  tempDict = dict()
  for rxn in metdict.keys():
    tempDict[rxn]=abs(metdict[rxn])

  #sort metdict by values
  import operator
  sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)

  metdict2 = dict()
  metdict2["Others-pos"]=0
  metdict2["Others-neg"]=0
  baseline = dict()
  pos_base=0
  neg_base=0
  i=0
  for TEMP in sorted_by_value:
    rxn = TEMP[0]
    if metdict[rxn]>0:
      if metdict[rxn] < total*threshold:
        if percentage:
          metdict2["Others-pos"]=metdict2["Others-pos"]+float(metdict[rxn]*100)/total
        else:
          metdict2["Others-pos"]=metdict2["Others-pos"]+metdict[rxn]
        continue
      base = pos_base
      if percentage:
        metdict2[rxn]=float(metdict[rxn]*100)/total
        pos_base = pos_base + float(metdict[rxn]*100)/total
      else:
        pos_base = pos_base + metdict[rxn]
        metdict2[rxn]=metdict[rxn]
    else:
      if abs(metdict[rxn]) < total*threshold:
        if percentage:
          metdict2["Others-neg"]=metdict2["Others-neg"]+float(metdict[rxn]*100)/total
        else:
          metdict2["Others-neg"]=metdict2["Others-neg"]+metdict[rxn]
        continue
      base = neg_base
      if percentage:
        metdict2[rxn]=float(metdict[rxn]*100)/total
        neg_base = neg_base + float(metdict[rxn]*100)/total
      else:
        neg_base = neg_base + metdict[rxn]
        metdict2[rxn]=metdict[rxn]
    i=i+1
    baseline[rxn]=base
  baseline["Others-pos"]=pos_base
  baseline["Others-neg"]=neg_base

  if show_plot:
    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.figure(figsize=(3,4))
    for rxn in metdict2.keys():
      if colourDict.keys().__contains__(rxn):
        plt.bar(1,metdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn])
      else:
        plt.bar(1,metdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
    plt.xlim(0.8,1.2)
    if percentage:
      plt.ylabel(metabolite+" produced/consumed (%)")
    else:
      plt.ylabel(metabolite+" produced/consumed (in moles)")
    handles, labels = plt.gca().get_legend_handles_labels()
    labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
    handles2=[handles[labels.index(i)] for i in labels2]
    lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')

  # create a pandas dataframe from the dictionary
  met_df = pd.DataFrame.from_dict(metdict2, orient='index', columns=['value'])

  return metdict2, met_df










#####

def estimateOutputFromNetCO2_Corinna_dev(model,netCO2uptake,Output_ID="Phloem_output_tx_light_dark",
                                      Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_light",
                                      CO2in_ID="CO2_tx_light",verbose=False, iterations=50, threshold=0.001, pFBA=True):
    

    
    """ 
    #SANU SHAMEER'S FUNCTION MODIFIED BY CORINNA HARTINGER
    
    this function compares the current predicted CO2 uptake rate to the user-defined CO2 uptake rate, and changes the bounds on the objective reaction until the 
    predicted value matches the user-defined value
    """

    temp_model = model.copy()
    
    # define objective reaction
    temp_model.objective = {temp_model.reactions.get_by_id(Output_ID): 1}
    
    
    #############################################################################################################################################################
    # test whether model can achieve desired CO2 uptake rate
    temp_model.reactions.get_by_id(CO2in_ID).bounds = (netCO2uptake, netCO2uptake)
    
    try:
        if pFBA:
            solution = cobra.flux_analysis.parsimonious.pfba(temp_model)
        else:
            solution = temp_model.optimize() #does not fail with error, just issues solver warning
            # would require some minimisation of difference expt. CO2 and achieved CO2, as is, the model doesn't decrease CO2 flux as objective flux is lowered
        print("Model can achieve desired CO2 rate, without discouraging CO2 refixation.")
        print(temp_model.reactions.get_by_id(Output_ID).flux)
    except:
        print("\033[91m"+"\nModel cannot achieve desired CO2 rate, even without discouraging CO2 refixation!\033[00m")
        # exit function
        raise Exception("\033[91m"+"\nModel cannot achieve desired CO2 rate, even without discouraging CO2 refixation!\033[00m")
    #############################################################################################################################################################
    
    
    # unconstrain upper bound of daytime CO2 uptake reaction
    temp_model.reactions.get_by_id(CO2in_ID).bounds = (-1000, 1000)
    
    #initial model optimisation
    try:
        if pFBA:
            solution = cobra.flux_analysis.parsimonious.pfba(temp_model)
        else:
            solution = temp_model.optimize()
    except:
        print("\033[91m"+"\nFirst optimisation: Model is infeasible!\033[00m")
        return
    else:
        print("First optimisation (unconstrained CO2 uptake): Model is feasible.")
        current_output = temp_model.reactions.get_by_id(Output_ID).flux

    print("Starting CO2 uptake =", temp_model.reactions.get_by_id(CO2in_ID).flux)
    print("Starting output flux:", Output_ID, "=", temp_model.reactions.get_by_id(Output_ID).flux)
    print("\n")
    
    
    # Create an empty dataframe to store the data
    import pandas as pd
    data = pd.DataFrame(columns=['Iteration', 'nCO2RC', 'CCE', 'CO2_in', 'CO2_out', 'Objective flux'])

    CCE = calculate_CCE(temp_model)
    nCO2RC = calculate_noct_CO2_refixation_coefficient(temp_model, verbose=False, threshold=0.00001, show_noflux=False)
    CO2_in = temp_model.reactions.CO2_tx_00.flux
    CO2_out = temp_model.reactions.CO2_tx_12.flux

    data = append_dict_to_df(data,{'Iteration': 0,  'nCO2RC': nCO2RC, 'CCE': CCE, 
                                'CO2_in': CO2_in, 'CO2_out': CO2_out,
                                'Objective flux': temp_model.reactions.get_by_id(Output_ID).flux})


    
    #Use a while loop to change Output flux until net CO2 rate is similar to given value 
    #(or loop counter hits limit)
    
    

    i=0
    while abs(temp_model.reactions.get_by_id(CO2in_ID).flux - netCO2uptake) > threshold and i<iterations:

        i += 1
        print("---- CO2 adjustment iteration " + str(i) + " ----")

        prev = temp_model.reactions.get_by_id(Output_ID).flux
        CCE = round((1 + (temp_model.reactions.CO2_tx_12.flux / temp_model.reactions.CO2_tx_00.flux)) * 100, 2)
        print("CCE:", CCE)

        # Obtain new bounds for Output flux by multiplying old bounds by 
        # ratio of user-defined C02 intake to current CO2 intake

        diff = temp_model.reactions.get_by_id(CO2in_ID).flux - netCO2uptake
        factor = (netCO2uptake / temp_model.reactions.get_by_id(CO2in_ID).flux)# ** 1.8 #add a power to the factor to make the change more extreme
        print("Factor:", round(factor, 4))

        now = prev * factor
        # now = prev + (prev*(netCO2uptake - temp_model.reactions.get_by_id(CO2in_ID).flux/netCO2uptake))

        print("prev:", prev)
        print("now:", now)

        temp_model.reactions.get_by_id(Output_ID).bounds = (now, now)

        try:
            if pFBA:
                solution = cobra.flux_analysis.parsimonious.pfba(temp_model)
            else:
                solution = temp_model.optimize()

            print("Optimisation successful.")

            current_output = now
            CCE = calculate_CCE(temp_model)
            nCO2RC = calculate_noct_CO2_refixation_coefficient(temp_model, verbose=False, threshold=0.00001, show_noflux=False)
            CO2_in = temp_model.reactions.CO2_tx_00.flux
            CO2_out = temp_model.reactions.CO2_tx_12.flux
        except:
            raise Exception("Iteration infeasible (output bounds =" + str(now) + ")")
            current_output = prev

            nCO2RC = 0
            CCE = 0
            CO2_in = 0
            CO2_out = 0


        finally:
            # Save the data to the dataframe
            data = append_dict_to_df(data,{'Iteration': i,  'nCO2RC': nCO2RC, 'CCE': CCE, 
                                'CO2_in': CO2_in, 'CO2_out': CO2_out,
                                'Objective flux': now})

    
        

        if verbose:
            print("Vc flux = "+str(temp_model.reactions.get_by_id(Vc_ID).flux))
            print("Current CO2 uptake = "+str(temp_model.reactions.get_by_id(CO2in_ID).flux))
            print("Target CO2 uptake = "+str(netCO2uptake))
            print("Output_ID before: "+str(prev))
            print("Output_ID after: "+str(now))
            CCE = round((1 +(temp_model.reactions.CO2_tx_12.flux/temp_model.reactions.CO2_tx_00.flux))*100,2)
            print("CCE:", CCE)
            print("Sum of fluxes:", solution.objective_value)
            print("\n")
    
    else:
        if verbose:
            print("Final results")
            print("----"+str(i)+"----")
            print("Vc flux = "+str(temp_model.reactions.get_by_id(Vc_ID).flux))
            print("Current CO2 uptake = "+str(temp_model.reactions.get_by_id(CO2in_ID).flux))
            print("Target CO2 uptake = "+str(netCO2uptake))
            print("Output_ID:", temp_model.reactions.get_by_id(Output_ID).flux)
            CCE = round((1 +(temp_model.reactions.CO2_tx_12.flux/temp_model.reactions.CO2_tx_00.flux))*100,2)
            print("CCE:", CCE)
            print("Sum of fluxes:", solution.objective_value)
            print("\n")
            
        if i < iterations:
            print("Net CO2 uptake value achieved.\n")
        else:
            raise Exception("\033[91m"+"\nMaximum iterations reached: Desired CO2 uptake value not reached after "+str(i)+" iterations.\033[00m")
        
        
    
    return solution, temp_model, data







def find_metabolite_rxns(model, metabolite_id):
    """
    This function takes a metabolic model and a metabolite ID as input and returns four lists of reaction IDs:
    1. A list of reaction IDs that involve the given metabolite.
    2. A list of reaction IDs that are bidirectional for the given metabolite.
    3. A list of reaction IDs that only consume the given metabolite.
    4. A list of reaction IDs that only produce the given metabolite.

    Parameters:
    model (cobra.Model): A metabolic model in COBRApy format.
    metabolite_id (str): Full or part of the ID of the metabolite to search for in the model.

    Returns:
    tuple: A tuple containing four lists of reaction IDs:
    1. A list of reaction IDs that involve the given metabolite.
    2. A list of reaction IDs that are bidirectional for the given metabolite.
    3. A list of reaction IDs that only consume the given metabolite.
    4. A list of reaction IDs that only produce the given metabolite.
    
    """
    met_rxns = []
    consume_rxns = []
    produce_rxns = []
    bidirectional_rxns = []
    
    for rxn in model.reactions:
        for met in rxn.metabolites:
            met_coefficient = rxn.get_coefficient(model.metabolites.get_by_id(met.id))
            
            if metabolite_id in met.id:
                met_rxns.append(rxn.id)
                
                if rxn.lower_bound < 0 and rxn.upper_bound > 0:
                    if rxn.id not in bidirectional_rxns:
                        bidirectional_rxns.append(rxn.id)
                
                elif (met_coefficient < 0 and rxn.lower_bound >= 0) or (met_coefficient > 0 and rxn.upper_bound <= 0):
                    if rxn.id not in consume_rxns:
                        consume_rxns.append(rxn.id)
                
                elif (met_coefficient < 0 and rxn.upper_bound <= 0) or (met_coefficient > 0 and rxn.lower_bound >= 0):
                    if rxn.id not in produce_rxns:
                        produce_rxns.append(rxn.id)
    
    print("#####  " +str(len(met_rxns)) +" reactions found to involve "+ metabolite_id+"  ####")
    print(str(len(bidirectional_rxns)) +" reactions found to be bidirectional for "+ metabolite_id)
    print(str(len(consume_rxns)) +" reactions found to only consume "+ metabolite_id)
    print(str(len(produce_rxns)) +" reactions found to only produce "+ metabolite_id)
    
    return (met_rxns, bidirectional_rxns, consume_rxns, produce_rxns)








def plot_fluxes_day_night(rxns, df_list):
    rxns_day = [r for r in rxns if r.endswith("_00")]
    rxns_night = [r for r in rxns if r.endswith("_12")]

    fig, ax = plt.subplots(nrows=len(rxns_day), ncols=1, figsize=(10, 30))

    colors = ["blue", "orange", "green", "red", "purple"]

    for count, (dayflux, nightflux) in enumerate(zip(rxns_day, rxns_night)):
        for i, df in enumerate(df_list):
            if sum(df.T[dayflux].dropna()) != 0:
                df.T.plot(x='CO2_tx_00', y=dayflux, kind='line', style="-.", ax=ax[count], label=f'df{i}', color=colors[i])
            if sum(df.T[nightflux].dropna()) != 0:
                df.T[nightflux] =(df.T[nightflux]*(-1))
                df.T.plot(x='CO2_tx_00', y=nightflux, kind='line', style="--", ax=ax[count], color=colors[i], label=f'df{i}')

        ax[count].set_xlabel('CO2_tx_00')
        ax[count].set_ylabel(dayflux.replace("_linker_00_to_12", ""))
        ax[count].set_title(dayflux.replace("_linker_00_to_12", "") +' vs CO2_tx_00')
        # ax[count].legend()

    plt.tight_layout()
    plt.show()





def plot_fluxes(rxns, df_list, figure_height=30):
    fig, ax = plt.subplots(nrows=len(rxns), ncols=1, figsize=(10, figure_height))

    colors = ["blue", "orange", "green", "red", "purple"]

    for count, flux in enumerate(rxns):
        for i, df in enumerate(df_list):
            if sum(df.T[flux].dropna()) != 0:
                df.T.plot(x='CO2_tx_00', y=flux, kind='line', style="-.", ax=ax[count], label=f'df{i}', color=colors[i])

        ax[count].set_xlabel('CO2_tx_00')
        ax[count].set_ylabel(flux.replace("_linker_00_to_12", ""))
        ax[count].set_title(flux.replace("_linker_00_to_12", "") +' vs CO2_tx_00')
        # ax[count].legend()

    plt.tight_layout()
    plt.show()



def calculate_sum_of_all_fluxes(solution):
    sum_of_fluxes = 0
    for item in solution:
        sum_of_fluxes += abs(item)

    return sum_of_fluxes

def photon_use(solution):
    solution["Photon_tx_00"]

    return photon_use





def plot_metabolite_budgets(df_name_list, all_models, solutions, metabolite, tag="_12"):
    import matplotlib.pyplot as plt
    import pandas as pd

    # create a DataFrame
    df_all = pd.DataFrame()
    dfs = []

    # create a list of data for the y-axis and add it to the DataFrame
    for i, model in enumerate(all_models):
        _, new_row_df = generateMetaboliteBudget(model, solutions[i], outfile="", show_plot=False, percentage=False,
                                                 day_or_night_tag=tag, save_plot_to="", threshold=0.02,
                                                 colourDict={},
                                                 metabolite=metabolite)
        new_row_df.columns = [df_name_list[i]]
        dfs.append(new_row_df)

    for i in range(len(dfs)):
        df_all = df_all.combine_first(dfs[i])
    colors = plt.cm.tab20(np.linspace(0, 1, len(df_all.index)))

    # create a stacked bar plot for each row in the DataFrame
    ax = df_all.T.plot(kind='bar', stacked=True, figsize=(10, 12), color=colors)

    # set the title and y-axis label for the plot
    ax.set_title(metabolite + ' produced/consumed by reaction')
    ax.set_ylabel(metabolite + ' produced/consumed')
    ax.axhline(y=0, color='black', linestyle='-')

    # display the plot
    plt.tight_layout()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

    return df_all


def plot_fluxes_obj_flux(rxns, df_list, objective_reaction = "AraCore_Biomass_tx_total", figure_height=30):
    fig, ax = plt.subplots(nrows=len(rxns), ncols=1, figsize=(10, figure_height))

    colors = ["blue", "orange", "green", "red", "purple", "brown", "turquoise", "pink", "olive", "cyan", "magenta", "yellow", "black"]

    for count, flux in enumerate(rxns):
        for i, df in enumerate(df_list):
            if sum(df.T[flux].dropna()) != 0:
                df.T.plot(x=objective_reaction, y=flux, kind='line', style="-.", ax=ax[count], label=f'df{i}', color=colors[i])

        ax[count].set_xlabel(objective_reaction)
        ax[count].set_ylabel(flux)
        ax[count].set_title(flux +' vs '+objective_reaction)
        # ax[count].legend()

    plt.tight_layout()
    plt.show()



def plot_flux_bars(df_name_list, all_models):
    """
    Plot flux bars for each model in all_models.

    Parameters:
    - df_name_list (list): List of labels for the x-axis.
    - all_models (list): List of models to plot.

    Returns:
    - None
    """

    # create a DataFrame
    df = pd.DataFrame()

    # create a list of data for the y-axis and add it to the DataFrame
    for i, model in enumerate(all_models):
        new_col = calc_carbon_flux(model, "AraCore_Biomass_tx")
        df.insert(i, df_name_list[i], new_col)

    # create a subplot for each row in the DataFrame
    fig, axs = plt.subplots(df.shape[0], figsize=(10, df.shape[0]*3))

    # plot bars for each column in the DataFrame
    for i, ax in enumerate(axs):
        bars = ax.bar(df.columns, df.iloc[i])
        ax.set_title(df.T.columns[i])
        ax.set_ylabel(df.T.columns[i])

        # set the color of the bars based on the label
        for bar, label in zip(bars, df.columns):
            if 'CETCH' in label:
                bar.set_color('orange')
            else:
                bar.set_color('lightblue')

    # display the plots
    plt.tight_layout()
    plt.show()



def find_metabolite_consumption_production_rxns(model, met_id, threshold = 0.0000001):
    metabolites = []

    for met in model.metabolites:
        if met_id in met.id:
            metabolites.append(met)

    met_rxns = {}

    for rxn in model.reactions:
        total_met = 0
        for met in rxn.metabolites:
            if met in metabolites:
                total_met += rxn.get_coefficient(met)
        if abs(total_met * rxn.flux) > threshold:
            met_rxns[rxn.id] = total_met * rxn.flux

    met_rxns_sorted = sorted(met_rxns.items(), key=lambda x: x[1])

    return met_rxns_sorted




def plot_trendline(df, x_column, y_column):
    # plot the data
    plt.plot(df.T[x_column], df.T[y_column], label=y_column)
    plt.xlabel(x_column)
    plt.ylabel(y_column)
    plt.title(f'{y_column} vs {x_column}')

    # calculate the slope
    x = df.T[x_column]
    y = df.T[y_column]
    x = x.dropna()
    y = y.dropna()
    slope, intercept = np.polyfit(x, y, 1)
    print(f'The slope of {y_column} vs {x_column} is:', slope)
    print(f'The intercept of {y_column} vs {x_column} is:', intercept)

    # plot the trendline
    plt.plot(x, slope*x + intercept, color='red', label='trendline', alpha=0.5)
    plt.legend()
    plt.show()

    return slope, intercept

def calc_met_turnover(model, mets, objective=""):
    """calculates the turnover of a metabolite in a model"""
    """Inputs: model (cobra model), mets (list of metabolite IDs)"""

    if objective != "":
        print("Results normalised by specfied reaction flux:", objective, model.reactions.get_by_id(objective).flux)
        print("\n")
    print("Turnover: \t Total \t Day \t Night")
    for met_id in mets:
        total_met = 0
        day_met = 0
        night_met= 0
        for rxn in model.reactions:
            temp_met_list = [met.id for met in rxn.metabolites]
            for item in temp_met_list:
                if item.startswith(met_id):
                    # print(rxn)
                    # print(rxn.flux)
                    # print("\n")
                    if rxn.flux * rxn.get_coefficient(item) > 0:
                        total_met += rxn.flux * rxn.get_coefficient(item)
                    if rxn.flux * rxn.get_coefficient(item) > 0 and rxn.id.endswith("_00"):
                        day_met += rxn.flux * rxn.get_coefficient(item)
                    elif rxn.flux * rxn.get_coefficient(item) > 0 and rxn.id.endswith("_12"):
                        night_met += rxn.flux * rxn.get_coefficient(item)
        
        if objective == "":
            print(met_id, "\t", total_met, "\t", day_met, "\t", night_met)
        else:
            total_met = total_met / model.reactions.get_by_id(objective).flux
            day_met = day_met / model.reactions.get_by_id(objective).flux
            night_met = night_met / model.reactions.get_by_id(objective).flux
            print(met_id, "\t", total_met, "\t", day_met, "\t", night_met)

    
def find_metabolite_producing_rxns(model, desired_metabolite):
    metabolite_producing_rxns = []
    metabolite_producing_rxns_dict = {"forward_rxns": [], "reverse_rxns": []}

    for rxn in model.reactions:
        metabolite_sum = 0
        fw = False
        rev = False
        for met in rxn.metabolites:
            for comp in model.compartments:
                if met.id == desired_metabolite + comp:
                    metabolite_sum += rxn.get_coefficient(met)
                    if (rxn.get_coefficient(met) > 0 and rxn.upper_bound > 0):
                        fw = True
                    elif (rxn.get_coefficient(met) < 0 and rxn.lower_bound < 0):
                        rev = True

        if metabolite_sum != 0:
            if fw:
                metabolite_producing_rxns_dict["forward_rxns"].append(rxn)
                metabolite_producing_rxns.append(rxn)
            elif rev:
                metabolite_producing_rxns_dict["reverse_rxns"].append(rxn)
                metabolite_producing_rxns.append(rxn)
        fw = False
        rev = False
    return metabolite_producing_rxns_dict



def set_reaction_constraint(model, desired_metabolite, objective_rxn):
    # create an empty reaction for total metabolite X-production

    temp_rxn = Reaction(objective_rxn, name=objective_rxn)
    temp_rxn.add_metabolites({})
    temp_rxn.bounds = (-1000, 1000)
    model.add_reactions([temp_rxn])
    
    coefficients = dict()

    rxns_dict = find_metabolite_producing_rxns(model, desired_metabolite)

    print("The following reactions were reversed in the model:")
    for rxn in rxns_dict["reverse_rxns"]:
        reaction = model.reactions.get_by_id(rxn.id)
        # print(rxn , rxn.bounds)
        for met in rxn.metabolites:
            met = model.metabolites.get_by_id(met.id)
            coeff = reaction.get_coefficient(model.metabolites.get_by_id(met.id))
            # print(coeff)
                                        
            rxn.add_metabolites({met: coeff * (-1)}) #to remove it
            rxn.add_metabolites({met: coeff * (-1)}) #to add it to the other side
        print(rxn.id , rxn.bounds)


    for fw_rxn in rxns_dict["forward_rxns"]:
        CO2_coeff = 0
        CO2_list = [CO2_met.id for CO2_met in fw_rxn.metabolites if CO2_met.id.startswith(desired_metabolite)]
        print(fw_rxn.id)
        print(CO2_list)
        for CO2_met in model.metabolites.query(desired_metabolite):
            if CO2_met.id in CO2_list:
                CO2_coeff = fw_rxn.get_coefficient(CO2_met.id)
        coefficients[fw_rxn.forward_variable] = 1 *CO2_coeff
        coefficients[fw_rxn.reverse_variable] = -1*CO2_coeff
        print(CO2_coeff)
    for rev_rxn in rxns_dict["reverse_rxns"]:
        CO2_coeff = 0
        CO2_list = [CO2_met.id for CO2_met in rev_rxn.metabolites if CO2_met.id.startswith(desired_metabolite)]
        print(rev_rxn.id)
        print(CO2_list)
        for CO2_met in model.metabolites.query(desired_metabolite):
            if CO2_met.id in CO2_list:
                CO2_coeff = rev_rxn.get_coefficient(CO2_met.id)

        coefficients[rev_rxn.forward_variable] = 1*CO2_coeff
        coefficients[rev_rxn.reverse_variable] = -1*CO2_coeff
        print(CO2_coeff)

    coefficients[model.reactions.get_by_id(objective_rxn).forward_variable] = -1
    coefficients[model.reactions.get_by_id(objective_rxn).reverse_variable] = 1

    constraint = model.problem.Constraint(0, lb=0, ub=0)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)

    return rxns_dict

    for fw_rxn in rxns_dict["forward_rxns"]:
        coefficients[fw_rxn.forward_variable] = 1
        # coefficients[fw_rxn.reverse_variable] = -1
    for rev_rxn in rxns_dict["reverse_rxns"]:
        # coefficients[rev_rxn.reverse_variable] = 1
        # coefficients[rev_rxn.forward_variable] = -1
        # coefficients[rev_rxn.reverse_variable] = -1
        coefficients[rev_rxn.forward_variable] = 1

    coefficients[model.reactions.get_by_id(objective_rxn).forward_variable] = -1
    coefficients[model.reactions.get_by_id(objective_rxn).reverse_variable] = 1

    constraint = model.problem.Constraint(0, lb=0, ub=0)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)

    return rxns_dict

def calculate_metabolite_production_consumption(model, metabolite_id):
    production = 0
    consumption = 0
    for rxn in model.reactions:
        if rxn.id[-2:] not in ("mc", "pc", "xc", "ec", "tx"):
            for met in rxn.metabolites:
                if met.id.startswith(metabolite_id):
                    coeff_flux = rxn.get_coefficient(model.metabolites.get_by_id(met.id)) * rxn.flux
                    if coeff_flux > 0:
                        production += abs(coeff_flux)
                        print("prod", rxn.id, round(rxn.flux, 5), round(coeff_flux, 5))
                    elif coeff_flux < 0:
                        consumption += abs(coeff_flux)
                        print("cons", rxn.id, round(rxn.flux, 5), round(coeff_flux, 5))
    print("prod:", round(production, 5))
    print("cons:", round(consumption, 5))

import pandas as pd

def create_dataframe(labels, solutions, models):
    # create a DataFrame
    df = pd.DataFrame()

    # create a list of data for the y-axis and add it to the DataFrame
    for i, model in enumerate(models):
        new_col = calc_carbon_flux(model,"AraCore_Biomass_tx")
        df.insert(i, labels[i], new_col)
 
    # create a DataFrame
    df3 = pd.DataFrame()

    # create a list of data for the y-axis and add it to the DataFrame
    for i, sol in enumerate(solutions):
        new_col = calculate_sum_of_all_fluxes(sol)
        df3.insert(i, labels[i], {"sum of fluxes":new_col})

    # create a DataFrame
    df2 = pd.DataFrame()

    # create a list of data for the y-axis and add it to the DataFrame
    for i, sol in enumerate(solutions):
        new_col = sol.fluxes["AraCore_Biomass_tx_total"]
        df2.insert(i, labels[i], {"objective flux total":new_col})

    # create a DataFrame
    df1 = pd.DataFrame()

    # create a list of data for the y-axis and add it to the DataFrame
    for i, sol in enumerate(solutions):
        new_col = sol.fluxes["Photon_tx_00"]
        df1.insert(i, labels[i], {"Photon flux":new_col})

    
    linker_list = ["MAL_v_linker_12_to_00", "CIT_v_linker_12_to_00", "STARCH_p_linker_00_to_12"]
    for linker in linker_list:
        df4 = pd.DataFrame()
        for i, sol in enumerate(solutions):
            
            new_col = sol.fluxes[linker]
            df4.insert(0, labels[i], {linker:new_col})
        df = df4.combine_first(df)

        # create a list of data for the y-axis and add it to the DataFrame

    # create a DataFrame
    df5 = pd.DataFrame()

    for i, sol in enumerate(solutions):
        new_col = sol.fluxes["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_00"]
        df5.insert(i, labels[i], {"Rubisco flux":new_col})

    df = df5.combine_first(df)
    df = df3.combine_first(df)
    df = df1.combine_first(df)
    df = df.combine_first(df2)
   

    return df


    ############################################################################################################################
# function to calculate reversiblity score according to Henry et al. 2009

def calc_rxn_reversibility_score(model, rxn_id):

    rxn = model.reactions.get_by_id(rxn_id)
    print(rxn)

    #check for low energy metabolites
    # sum up number of low energy metabolites * number of its molecules
    low_E_mets = ["CARBON_DIOXIDE", "HCO3", "CO_A", "Pi", "PPI", "ACP"]
    low_E_score = 0
    for met in rxn.metabolites:
        for item in low_E_mets:
            if met.id.startswith(item):
                # n = abs(rxn.get_coefficient(met.id))
                n = rxn.get_coefficient(met.id)
                low_E_score +=  1 * n
    print(low_E_score)

    #check for ATP/ADP/AMP

    min_ATP_ADP_Pi = []
    min_ATP_AMP_PPI = []

    for met in rxn.metabolites:
        if met.id.startswith("ATP"):
            min_ATP_ADP_Pi.append(abs(rxn.get_coefficient(met.id)))
            min_ATP_AMP_PPI.append(abs(rxn.get_coefficient(met.id)))
        elif met.id.startswith("ADP"):
            min_ATP_ADP_Pi.append(abs(rxn.get_coefficient(met.id)))
        elif met.id.startswith("Pi"):
            min_ATP_ADP_Pi.append(abs(rxn.get_coefficient(met.id)))
        elif met.id.startswith("AMP"):
            min_ATP_AMP_PPI.append(abs(rxn.get_coefficient(met.id)))
        elif met.id.startswith("PPI"):
            min_ATP_AMP_PPI.append(abs(rxn.get_coefficient(met.id)))

    if len(min_ATP_ADP_Pi) == 0:
        min_ATP_ADP_Pi = [0,0,0] # if non of ATP/ADP/Pi was present, initiate a list with 0s
    elif len(min_ATP_ADP_Pi) != 0 and len(min_ATP_ADP_Pi) != 3:
        min_ATP_ADP_Pi.append(0)

    if len(min_ATP_AMP_PPI) == 0:
        min_ATP_AMP_PPI = [0,0,0]
    elif len(min_ATP_AMP_PPI) != 0 and len(min_ATP_AMP_PPI) != 3: #if the list is not empty, but also not full, append 0 because one of the above metabolites is not present. otherwise the minimum is not defined correctly
        min_ATP_AMP_PPI.append(0)

    print(min_ATP_ADP_Pi)
    print(min_ATP_AMP_PPI)

    #calculate reversiblity score
    rev_score = min(min_ATP_ADP_Pi) + min(min_ATP_AMP_PPI) - low_E_score
    print("Srev =", str(min(min_ATP_ADP_Pi)),"+", str(min(min_ATP_AMP_PPI)), "-",str(low_E_score))
        
    return rev_score
############################################################################################################################



### function to plot the day- and nighttime energy and reducing power demand of added pathway

def plot_day_night_pathway_energy_demand(whole_model, pathway_model, day_tag="_00", night_tag="_12"):
    """determines the names of a the orthologous pathway from the pathway_model.reactions ids;
    finds the reactions of the pathway inside the whole_model, with both the day_tag and the night_tag;
    calculates the ATP consumption from these reactions;
    plots the ATP consumption for day_tag pathway reactions and night_tag pathways reactions separated on the same plot"""
    
    # Get the orthologous pathway reactions from the pathway_model
    pathway_reactions = pathway_model.reactions
    
    # Find the corresponding reactions in the whole_model with the day_tag and night_tag
    day_reactions = []
    night_reactions = []
    for reaction in pathway_reactions:
        for rxn in whole_model.reactions:
            if reaction.id in rxn.id and day_tag in rxn.id:
                day_reactions.append(whole_model.reactions.get_by_id(rxn.id))
            elif reaction.id in rxn.id and night_tag in rxn.id:
                night_reactions.append(whole_model.reactions.get_by_id(rxn.id))
        

    
    # Calculate the ATP consumption for day_tag pathway reactions
    # day_atp_consumption = sum(abs(reaction.flux) for reaction in day_reactions if reaction.flux < 0)
    day_atp_consumption = 0
    for rxn in day_reactions:
        for met in rxn.metabolites:
            if "ATP_" in met.id:
                day_atp_consumption += rxn.get_coefficient(met.id) * rxn.flux

    # Calculate the ATP consumption for night_tag pathway reactions
    # night_atp_consumption = sum(abs(reaction.flux) for reaction in night_reactions if reaction.flux < 0)
    night_atp_consumption = 0
    for rxn in night_reactions:
        for met in rxn.metabolites:
            if "ATP_" in met.id:
                night_atp_consumption += rxn.get_coefficient(met.id) * rxn.flux
                print(rxn)
                print(rxn.flux)
                continue

    print(day_atp_consumption, night_atp_consumption)
    
    # Plot the ATP consumption for day_tag and night_tag pathway reactions
    import matplotlib.pyplot as plt
    plt.bar(["Day", "Night"], [day_atp_consumption, night_atp_consumption])
    plt.xlabel("Time of Day")
    plt.ylabel("ATP Consumption")
    plt.title("Day-Night Pathway Energy Demand")
    plt.show()



###########
import networkx as nx
from pyvis.network import Network

def visualize_model(model, exclude_met):
    # Create an empty graph
    G = nx.Graph()

    # Add nodes for metabolites
    for metabolite in model.metabolites:
        if metabolite.id not in exclude_met:
            G.add_node(metabolite.id, type='metabolite', color='orange')

    # Add nodes for reactions
    for reaction in model.reactions:
        G.add_node(reaction.id, type='reaction', color='blue')

    # Create edges depending on the metabolites in each reaction in the model
    for reaction in model.reactions:
        for metabolite in reaction.metabolites:
            if metabolite.id not in exclude_met:
                G.add_edge(reaction.id, metabolite.id)

    # Visualize the graph as a force-directed graph
    pos = nx.spring_layout(G)
    node_colors = ['orange' if G.nodes[node]['type'] == 'metabolite' else 'lightblue' for node in G.nodes]
    nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=500, font_size=8)
    plt.show()

    # Create an interactive graph using pyvis
    nt = Network(notebook=True)
    nt.from_nx(G)
    nt.save_graph(str(model.name)+".html")



#################################################################################################################################
def extract_mini_network(seed_rxns, reaction_degrees, model, exclude_met):
    seed_rxn_dict = {}
    for seed_rxn in seed_rxns:

        current_rxns = [seed_rxn]
        next_rxns = []
        mini_network_rxns = []

        for i in range(reaction_degrees):
            for reaction in current_rxns:
                if reaction in mini_network_rxns:
                    continue

                rxn = model.reactions.get_by_id(reaction)

                for met in rxn.metabolites:
                    met = model.metabolites.get_by_id(met.id)
                    if met.id not in exclude_met:
                        for rxn2 in met.reactions:
                            if rxn2.id not in mini_network_rxns and rxn2.id not in next_rxns and rxn2.id not in current_rxns:
                                next_rxns.append(rxn2.id)

            [mini_network_rxns.append(rxn) for rxn in current_rxns]
            print("-----", i, "degrees", "-----")
            print(len(current_rxns))
            current_rxns = next_rxns
            next_rxns = []

        
        print("total number of reactions after searching for", i, "degrees of reactions:", len(mini_network_rxns))
        print("\n")
        seed_rxn_dict[seed_rxn] = len(mini_network_rxns)

    return seed_rxn_dict




################ copied from https://stackoverflow.com/questions/70837397/good-alternative-to-pandas-append-method-now-that-it-is-being-deprecated
def append_dict_to_df(df, dict_to_append):
    df = pd.concat([df, pd.DataFrame.from_records([dict_to_append])])
    return df

    """
    use:

    # Creating an empty dataframe
    df = pd.DataFrame(columns=['a', 'b'])

    # Appending a row
    df = append_dict_to_df(df,{ 'a': 1, 'b': 2 })
    """


#################################################################################################################################
def plot_eO_iteration_results(df):
    """
    Plots the objective flux, CO2_in, and nCO2RC against iterations.
    To be used after the estimateOutput function.

    Parameters:
    df (pandas.DataFrame): The dataframe containing the data to be plotted.

    Returns:
    None
    """
    fig, ax1 = plt.subplots()

    # Plot objective flux on the left y-axis with decreased opacity and square shape
    ax1.scatter(df['Iteration'], df['Objective flux'], color='blue', alpha=0.5, marker='s')
    ax1.plot(df['Iteration'], df['Objective flux'], color='blue', alpha=0.3)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Objective Flux', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.set_ylim(0, 0.135)

    # Create a twin Axes object for CO2 on the right y-axis
    ax2 = ax1.twinx()
    ax2.set_ylim(0, 12.5)

    # Plot CO2_in on the right y-axis with decreased opacity and circle shape
    ax2.scatter(df['Iteration'], df['CO2_in'], color='green', alpha=0.5, marker='o')
    ax2.plot(df['Iteration'], df['CO2_in'], color='green', alpha=0.3)
    ax2.set_ylabel('CO2 flux', color='black')

    ax2.tick_params(axis='y', labelcolor='green')
    # Add a dotted line at ax2 = 6
    ax2.axhline(y=6.643, color='gray', linestyle='dotted', label='expt. CO2 uptake')

    ## Plot CO2_out on the right y-axis with decreased opacity and triangle shape
    #ax2.scatter(df['Iteration'], -df['CO2_out'], color='magenta', alpha=0.5, marker='^')
    #ax2.plot(df['Iteration'], -df['CO2_out'], color='magenta', alpha=0.3)

    # Create a third Axes object for nCO2RC on the right y-axis
    #ax3 = ax1.twinx()
    #ax3.spines['right'].set_position(('outward', 60))  # Move the third y-axis to the right
    #ax3.set_ylim(0, 1.5)

    # Plot nCO2RC on the right y-axis with decreased opacity and diamond shape
    #ax3.scatter(df['Iteration'], df['nCO2RC'], color='magenta', alpha=0.5, marker='D')
    #ax3.plot(df['Iteration'], df['nCO2RC'], color='magenta', alpha=0.3)
    #ax3.set_ylabel('nCO2RC', color='magenta')

    #ax3.tick_params(axis='y', labelcolor='magenta')
    #ax3.set_ylim(0, 1)

    plt.title('Objective Flux and CO2_in vs Iterations')
    plt.show()


##############################
def plot_night_carbon_fluxes(metrics_df):
    subset_df = metrics_df.T[['night carbon biomass flux', 'CO2 nightime exchange', 'Carbon night to day', '% nCO2R']]

    # Define custom colors for each item
    colors = ['green', 'red', 'blue', 'orange']

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    subset_df[['night carbon biomass flux', 'CO2 nightime exchange', 'Carbon night to day']].plot.bar(stacked=True, ax=ax1, color=colors[:3])
    subset_df['% nCO2R'].plot(style='o', ax=ax2, color=colors[3])  # Plot nCO2RC as a dot
    ax2.set_ylim(0, 105)

    ax1.set_ylabel('Carbon flux (umol/m2/s)')
    ax2.set_ylabel('noct. CO2 refixation (%)')

    ax1.set_title('Nocturnal carbon fluxes')

    ax1.hlines(0, -1, 4, colors='black')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    plt.show()


##################################
    
def save_model_and_fluxes_json(model, solution, file_name):
    import json

    # save the model as JSON
    cobra.io.save_json_model(model, file_name+".json")

    # save a single flux vector as JSON
    with open(file_name+"_fluxes.json", 'w') as f:
        json.dump(solution.fluxes.to_dict(), f)