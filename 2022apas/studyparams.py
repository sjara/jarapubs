"""
This file contains the default names and queries used for the study on
Passive exposure effects on learning (active vs active+passive)
"""

STUDY_NAME = '2022apas'


# -- List of animals used in this study --

# Excluding 'pamo020' who took a while to pass criteria 2->3
MICE_ACTIVE_ONLY_COH2 = ['pamo009', 'pamo011', 'pamo012', 'pamo013', 'pamo016',
                         'pamo019', 'pamo022', 'pamo026'] # 'pamo020'
MICE_ACTIVE_PASSIVE = ['pamo010', 'pamo014', 'pamo015', 'pamo017', 'pamo018',
                       'pamo021', 'pamo023', 'pamo025', 'pamo024']

MICE_ALL_COH2 = MICE_ACTIVE_ONLY_COH2 + MICE_ACTIVE_PASSIVE
MICE_ALL_COH2.sort()

MICE_ACTIVE_ONLY_COH3 = ['pamo027', 'pamo029','pamo031', 'pamo034', 'pamo037', 'pamo038']
MICE_PASSIVE_THEN_ACTIVE = ['pamo028', 'pamo030', 'pamo032', 'pamo033', 'pamo035',
                            'pamo036', 'pamo039', 'pamo040', 'pamo041']
MICE_ALL_COH3 = ['pamo027', 'pamo028', 'pamo029', 'pamo030', 'pamo031', 'pamo032',
                 'pamo033', 'pamo034', 'pamo035', 'pamo036', 'pamo037', 'pamo038',
                 'pamo039', 'pamo040', 'pamo041']

MICE_ALL_COH4 = ['pamo042', 'pamo043', 'pamo044', 'pamo045', 'pamo046', 'pamo047',
                 'pamo048', 'pamo049', 'pamo050', 'pamo051', 'pamo052', 'pamo053',
                 'pamo054', 'pamo055', 'pamo056', 'pamo057', 'pamo058', 'pamo059']
MICE_ACTIVE_ONLY_COH4 = MICE_ALL_COH4


MICE_ACTIVE_ONLY_ALL = MICE_ACTIVE_ONLY_COH2 + MICE_ACTIVE_ONLY_COH3 + MICE_ACTIVE_ONLY_COH4
#MICE_ACTIVE_ONLY_ALL = MICE_ACTIVE_ONLY_COH2 + MICE_ACTIVE_ONLY_COH4

SESSIONS_RANGE_COH2 = ['2022-01-14', '2022-04-02']  # Cohort 2, all stages
SESSIONS_RANGE_COH3 = ['2022-04-22', '2022-06-26']  # Cohort 3, all stages 
SESSIONS_RANGE_COH4 = ['2022-06-28', '2022-08-14']  # Cohort 4, all stages 


SESSIONS_COH2_STAGE2 = ['2022-01-20', '2022-01-31'] # Cohort 2, Stage 2
SESSIONS_COH3_STAGE2 = ['2022-04-28', '2022-05-06'] # Cohort 3, Stage 2
SESSIONS_COH4_STAGE2 = ['2022-07-04', '2022-07-12'] # Cohort 3, Stage 2


# -- Exclude last session stage 3 (failed rigs for COH2) and extra session for very slow mouse --
SESSIONS_COH2_STAGE3 = ['2022-02-01', '2022-02-26'] # Cohort 2, Stage 3
SESSIONS_COH3_STAGE3 = ['2022-05-07', '2022-06-01'] # Cohort 3, Stage 3
SESSIONS_COH4_STAGE3 = ['2022-07-13', '2022-08-07'] # Cohort 3, Stage 3


LAST_SESSIONS_COH2_STAGE3 = ['2022-02-23', '2022-02-26'] # Last 4 sessions, Cohort 2, Stage 3
LAST_SESSIONS_COH3_STAGE3 = ['2022-05-28', '2022-05-31'] # Last 4 sessions, Cohort 3, Stage 3
LAST_SESSIONS_COH4_STAGE3 = ['2022-08-04', '2022-08-07'] # Last 4 sessions, Cohort 3, Stage 3

#EARLY_SESSIONS_COH2_STAGE4 = ['2022-02-28', '2022-03-03'] # First 4 sessions, Cohort 2, Stage 4
#LATE_SESSIONS_COH2_STAGE4 = ['2022-03-14', '2022-03-17']  # Day 15-18, Cohort 2, Stage 4
#VLATE_SESSIONS_COH2_STAGE4 = ['2022-03-20', '2022-03-30']  # Day 21-31, Cohort 2, Stage 4

# -- Check the wiki page for this study to see what happened on each session --
EXCLUDE_SESSIONS = {'pamo022':['2022-03-13'],
                    'pamo023':['2022-01-24'],
                    'pamo026':['2022-01-21'],
                    'pamo032':['2022-05-05'],
                    'pamo039':['2022-05-13'],
                    'pamo052':['2022-07-03'],
                    'pamo048':['2022-07-07'],
                    'pamo051':['2022-07-10'],
                    'pamo056':['2022-07-23']}

# -- Sessions with unexplicably bad behavior (probably rig issues?) --
UNRELIABLE_SESSIONS = {'pamo011':['2022-02-21'],
                       'pamo014':['2022-02-11'],
                       'pamo015':['2022-02-21'],
                       'pamo024':['2022-02-14'],
                       'pamo037':['2022-05-28']}
