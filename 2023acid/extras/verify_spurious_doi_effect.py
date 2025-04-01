"""
Verify that the effecst of DOI on oddball enhancement are not
simply the result of an overall decrease in response.
"""

# run figure_oddball_enhancement.py

stim

selcells = (responsive[stim] & steadyOEI & enoughTrials)

# These would be for DOI (the last 'reagent')
stnd = standardEvokedFR[selcells]
odd = oddballEvokedFR[selcells]

median(stnd)
median(odd)

origOEI = studyutils.modulation_index(odd, stnd)

# Evaluate what happens to the OEI if all firing rates are reduced by 50%
factor = 0.5
stndNew = stnd * factor
oddNew = odd * factor

newOEI = studyutils.modulation_index(oddNew, stndNew)

print(f'OEI: {median(origOEI):0.4f} -> {median(newOEI):0.4f} (factor={factor})')

