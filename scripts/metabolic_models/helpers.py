
from cobra import Reaction


def convertReactionsUnidirectional(cobra_model):
    new_reactions = []

    # Iterate through each reaction in the model
    for reaction in cobra_model.reactions:
        if reaction.lower_bound < 0:
            # Create a new reaction
            new_reaction = Reaction(reaction.id + "_rev", name = reaction.name + "_rev")  # Name the new reaction
            new_reaction.lower_bound = 0  # Set appropriate bounds (you might want to adjust this)
            new_reaction.upper_bound = -reaction.lower_bound  # You may choose to set it to a similar upper bound
            new_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            
            # Switch the reactants and products
            new_reaction.add_metabolites({metabolite: -coef for metabolite, coef in zip(reaction.reactants, [reaction.metabolites[i] for i in reaction.reactants])})
            new_reaction.add_metabolites({metabolite: -coef for metabolite, coef in zip(reaction.products, [reaction.metabolites[i] for i in reaction.products])})

            # Add the new reaction to the model
            cobra_model.reactions.get_by_id(reaction.id).lower_bound = 0
            new_reactions.append(new_reaction)
    cobra_model.add_reactions(new_reactions)

    return(cobra_model)