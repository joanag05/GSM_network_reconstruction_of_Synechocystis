from cobra import Model


def get_reactions(id_lipid, reactions_list, visited):
    lipid = model_cyano.metabolites.get_by_id(id_lipid + "__cytop")
    for reaction in lipid.reactions:
        if lipid in reaction.products or reaction.lower_bound < 0:
            reactions_list.add(reaction)
            for reactant in reaction.reactants:
                if reactant.id not in visited and (reactant.id.startswith("BMGC") or reactant.id.startswith("mgdg") or reactant.id.startswith("dgdg") or reactant.id.startswith("sqdg")):
                    visited.add(reactant.id)
                    reactions_to_add = get_reactions(reactant.id.split("__")[0], reactions_list, visited)
                    for reaction_to_add in reactions_to_add:
                        reactions_list.add(reaction_to_add)
    return reactions_list



from cobra.io import read_sbml_model, write_sbml_model

model_cyano = read_sbml_model('model_lipids_cyano.xml')


id_lipid = ["mgdg140183", "mgdg160184", "mgdg160183", "mgdg160182", "mgdg160181", "mgdg182183", "mgdg181183", "mgdg160200", "mgdg190210",
            ]

id_metabolites = [metabolite.id for metabolite in model_cyano.metabolites]
results = set()

for lipid in id_lipid:
    if lipid + "__cytop" in id_metabolites:
        reactions = get_reactions(lipid, set(), set())
        for reaction in reactions:
            for met in reaction.metabolites:
                met.compartment = "cytop"
            results.add(reaction)

new_model = Model()

new_model.add_reactions(results)

# write
write_sbml_model(new_model, "model_lipids_cyano_yo.xml")
