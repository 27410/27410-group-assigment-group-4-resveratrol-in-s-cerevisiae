def read_yeast_model(model_path):
    """This function is a modified version of the one available in the 
    GitHub repository yest-GEM (https://github.com/SysBioChalmers/yeast-GEM).
    It was developed by researchers in the department of Systems and Synthetic Biology 
    at Chalmers University og Techonology.

    Lu, H. et al. A consensus S. cerevisiae metabolic model Yeast8 and its ecosystem for 
    comprehensively probing cellular metabolism. Nature Communications 10, 3586 (2019). 
    doi:10.1038/s41467-019-11581-3

    Moreover, this function has been preivously modified in the year 2021 by students at
    DTU. Their code was use exclusively in this section. 
    See group-assignment-2021-group_10_rhb_cerevisiae.

    The function imports an SBML model and after imposing the necessary changes using COBRA, 
    exports it as BiGG complaint. 

    ----------
    make_bigg_compliant : bool, optional
        Whether the model should be initialized with BiGG compliance or not.
        If false, the original ids/names/compartments will be used instead.
    Returns
    -------
    cobra.core.Model
    """
    make_bigg_compliant=True

    import csv
    from cobra.io import read_sbml_model, write_sbml_model
    from copy import copy
    from dotenv import find_dotenv
    import os
    from os.path import dirname
    
    # Load model:
    MODEL_PATH = model_path
    model = read_sbml_model(MODEL_PATH)
    
    dotenv_path = find_dotenv()
    REPO_PATH = dirname(dotenv_path)

    # Check if already BiGG compliant:
    is_bigg_compliant = "x" in model.compartments

    # Convert to BiGG compliant if not already:
    if not is_bigg_compliant and make_bigg_compliant:
        # Load met/rxn dictionaries:
        def load_bigg_dict(bigg_file_path):
            bigg_dict = {}
            with open(bigg_file_path) as bigg_file:
                bigg_reader = csv.reader(bigg_file, delimiter=',')
                for row in bigg_reader:
                    bigg_dict[row[0]] = row[1]
            return bigg_dict
        data_path = f"{REPO_PATH}/data/databases"
        met_bigg_dict = load_bigg_dict(f"data/databases/BiGGmetDictionary_newIDs.csv")
        rxn_bigg_dict = load_bigg_dict(f"data/databases/BiGGrxnDictionary_newIDs.csv")

        # Function for adding unique ids to the model:
        def add_new_id(model_element, new_id):
            original_id = copy(new_id)
            id_assigned = False
            copy_number = 1
            while not id_assigned:
                try:
                    if hasattr(model_element, "compartment"):  # metabolites
                        model_element.id = f"{new_id}_{model_element.compartment}"
                    else:  # reactions
                        model_element.id = new_id
                    id_assigned = True
                except ValueError:
                    new_id = f"{original_id}_copy{str(copy_number)}"
                    copy_number += 1

        # Metabolite changes:
        comp_dic = {"er":"r", "erm":"rm", "p":"x"}
        for met in model.metabolites:
            # Save original id in notes:
            met.notes["Original ID"] = met.id
            # Change name to not include compartment at the end:
            met.name = met.name.replace(f" [{model.compartments[met.compartment]}]", "")
            # Change compartment info:
            if met.compartment in comp_dic:
                met.compartment = comp_dic[met.compartment]
            # Update id with BiGG information:
            if "bigg.metabolite" in met.annotation:
                add_new_id(met, met.annotation['bigg.metabolite'])
            elif met.id in met_bigg_dict:
                add_new_id(met, met_bigg_dict[met.id])
            else:
                met.id = met.id.replace(f"[{met.compartment}]", f"_{met.compartment}")

        # Compartment changes:
        comps = model.compartments
        comps["r"] = "endoplasmic reticulum"
        comps["rm"] = "endoplasmic reticulum membrane"
        comps["x"] = "peroxisome"
        model.compartments = comps

        # Reaction changes:
        for rxn in model.reactions:
            # Update id with BiGG information:
            if "bigg.reaction" in rxn.annotation:
                rxn.notes["Original ID"] = rxn.id
                add_new_id(rxn, rxn.annotation['bigg.reaction'])
            elif rxn.id in rxn_bigg_dict:
                rxn.notes["Original ID"] = rxn.id
                add_new_id(rxn, rxn_bigg_dict[rxn.id])

    return model


def write_yeast_model(model, model_output):
    from cobra.io import read_sbml_model, write_sbml_model
    """Writes the SBML file of the yeast model using COBRA.
    Parameters
    ----------
    model : cobra.core.Model
        Yeast model to be written
    model_output: name of the model
    """
    write_sbml_model(model, model_output)