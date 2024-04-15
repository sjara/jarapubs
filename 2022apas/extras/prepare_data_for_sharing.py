"""
A total of 27 wild-type C57BL/6J adult mice were used in this study.

Create a zip file of code with:
rm -r /var/tmp/Schmid2023/code/__pycache__
cd /var/tmp/
zip -r Schmid2023.zip Schmid2023 -x "*~"
"""

import sys
sys.path.append('..')
import studyparams
import studyutils

subjects = (studyparams.MICE_ACTIVE_ONLY_COH2 + ['pamo020'] +
            studyparams.MICE_ACTIVE_PASSIVE +
            studyparams.MICE_PASSIVE_THEN_ACTIVE)

# -- Copy relevant behavior data --
if 0:
    for subject in subjects:
        print(f'rsync -av /data/behavior/{subject} /var/tmp/Schmid2023/behavior_data/')
# You should exclude passive sessions from pamo010 and pamo028 (*/pamo*p?.h5')
# and a few wrong sessions (files have longer names)

if 0:        
    allmicestr = ("mice = {" +
                  "'active_only': ['" + "', '".join(studyparams.MICE_ACTIVE_ONLY_COH2) + "'], " +
                  "'unsuccessful': ['pamo020'], " +
                  "'active_plus_passive': ['" + "', '".join(studyparams.MICE_ACTIVE_PASSIVE) + "'], " +
                  "'passive_then_active': ['" + "', '".join(studyparams.MICE_PASSIVE_THEN_ACTIVE) + "']} " )
    print(allmicestr)

if 0:
    mice = {'active_only': ['pamo009', 'pamo011', 'pamo012', 'pamo013',
                            'pamo016', 'pamo019', 'pamo022', 'pamo026'],
            'unsuccessful': ['pamo020'],
            'active_plus_passive': ['pamo010', 'pamo014', 'pamo015', 'pamo017', 'pamo018',
                                    'pamo021', 'pamo023', 'pamo025', 'pamo024'],
            'passive_then_active': ['pamo028', 'pamo030', 'pamo032', 'pamo033', 'pamo035',
                                    'pamo036', 'pamo039', 'pamo040', 'pamo041']} 


#psyCurveData = studyutils.load_stage4('late')

if 0:
    stages = [0,1,2,3,4]
    sessions = {}

    for subject in subjects:
        sessions[subject] = {}
        for stage in stages:
            sessionThisStage = studyutils.get_sessions(subject, stage=stage)
            sessions[subject][stage] = sessionThisStage
    np.savez('/tmp/sessions_each_stage.npz', sessions=sessions)

    
