"""
Number of passive trials.

Because the data is only saved for one animal (one name), I'm instead relying on
data from the spreadsheet by Haziq.

pamo animals passive exposure session records and other related documents
https://docs.google.com/spreadsheets/d/10Is9ZeAYlhHc2juNxRGxmKeKBhFXvYt-DUrJCbJcfDQ/edit#gid=1438287531

"""

totalCoh2 = 192866
nSessionsCoh2 = 55 # Excluding the day they got 0

totalCoh3 = 81101
nSessionsCoh3 = 22

trialsPerSessionCoh2 = totalCoh2/nSessionsCoh2
trialsPerSessionCoh3 = totalCoh3/nSessionsCoh3

trialsPerSession = (totalCoh2+totalCoh3)/(nSessionsCoh2+nSessionsCoh3)

print(f'trialsPerSessionCoh2: {trialsPerSessionCoh2}')
print(f'trialsPerSessionCoh3: {trialsPerSessionCoh3}')
print(f'trialsPerSession: {trialsPerSession}')

