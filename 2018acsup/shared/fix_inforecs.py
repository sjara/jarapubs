"""
Update inforecs to new versions (post 2021).
"""

import os

command1 = """sed -i "s/, info/, probe='A4x2-tet', info/" {}"""
command2 = """sed -i "s/tetrodes/egroups/" {}"""

inforecPath = '/home/sjara/src/jarainfo/inforecordings/'

MICE = ['band016', 'band027', 'band028', 'band029', 'band030', 'band031', 'band034','band037','band038','band044','band045','band054','band059','band060']
#MICE = ['band005', 'band015']

for subject in MICE:
    print(f'Processing {subject}')
    filename = os.path.join(inforecPath, f'{subject}_inforec.py')
    if 0:
        os.system(command1.format(filename)) # Run this only once, ever!
    os.system(command2.format(filename))
    #print(command1.format(filename))
   

