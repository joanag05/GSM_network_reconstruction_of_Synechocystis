from cobra import Model
from cobra.manipulation import remove_genes


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
            "dgdg161183",
            "dgdg160183",
            "dgdg160182",
            "dgdg182183",
            "sqdg141160",
            "sqdg140162",
            "sqdg140183",
            "sqdg160162",
            "sqdg160161",
            "sqdg160160",
            "sqdg160172",
            "sqdg160171",
            "sqdg160170",
            "sqdg160184",
            "sqdg160183",
            "sqdg160182",
            "sqdg160181",
            "sqdg160180",
            "sqdg170183",
            "sqdg170182",
            "sqdg160190",
            "sqdg160204",
            "sqdg160203",
            "sqdg180182",
            "sqdg180181",
            "sqdg160200",
            "BMGC427",
            "BMGC682419",
            "BMGC14247",
            "BMGC428",
            "BMGC316",
            "BMGC315",
            "BMGC7294",
            "BMGC14561",
            "BMGC7455",
            "BMGC14238",
            "BMGC14239",
            "BMGC86142",
            "BMGC99059",
            "BMGC71918", "BMGC5629", "BMGC80156", "BMGC113284", "BMGC98263", "BMGC99055", "BMGC74677", "BMGC106626",
            "BMGC101001"
            ]

id_metabolites = [metabolite.id for metabolite in model_cyano.metabolites]
results = set()
not_in_model = set()

for lipid in id_lipid:
    if lipid + "__cytop" in id_metabolites:
        reactions = get_reactions(lipid, set(), set())
        for reaction in reactions:
            for met in reaction.metabolites:
                met.compartment = "cytop"
            results.add(reaction)
    else:
        print(lipid + "__cytop" + " not in model")
        not_in_model.add(lipid)

with open("not_in_model.txt", "w") as f:
    for lipid in not_in_model:
        f.write(lipid + "\n")

new_model = Model()

new_model.add_reactions(results)

groups_to_add = []
for group in model_cyano.groups:
    if any([reaction in results for reaction in group.members]):
        to_remove = []
        for reaction in group.members:
            if reaction not in results:
                to_remove.append(reaction)
        group.remove_members(to_remove)
        groups_to_add.append(group)

new_model.add_groups(groups_to_add)
# write

remove_genes(new_model, new_model.genes, remove_reactions=False)

write_sbml_model(new_model, "model_lipids_cyano_yo.xml")
