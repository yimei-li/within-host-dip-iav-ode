# Data

CSV behind the figures. Columns are in each file's header row.
Figure and table numbers follow the current manuscript.

invivo_observations.csv: mouse V/D/IFN. Fig 2, 4, S6, S7, S8.
invitro_observations.csv: A549 V/D/IFN. Fig 1D, S1, S2.
reporter_doseresponse.csv: reporter assay. Fig 1A-C.
pony_observations.csv: pony virus/IFN. Fig S4.
baccam_patient1_observations.csv: patient-1 titer. Fig S5.
posteriors_invivo_withDIP.csv: with-DIP fit medians. Table S1; Fig 4, 5, S3, S8, S9.
posteriors_invivo_withDIP_noFeedback.csv: feedback-off medians. Fig 4, feedback-off branch.
posteriors_invivo_noDIP.csv: no-DIP fit medians. Table S2; Fig S3, S8.
posteriors_invitro_withDIP.csv: in vitro with-DIP. Fig S2.
posteriors_invitro_noDIP.csv: in vitro no-DIP. Fig S1.
posteriors_joint_noDIP_sharedcumF50.csv: joint, shared cumF50. Fig S6.
posteriors_joint_noDIP_groupcumF50.csv: joint, per-group cumF50. Fig S7.
pony_map_params.csv: pony fit values. Fig S4.
baccam_fit_params.csv: patient-1 fit values. Fig S5.
constants.csv: fixed inputs (T0, m1, F0, LOD, nu_F).

The pairwise posterior figures in ../figures/ come from the raw draws of the same
with-DIP in vivo fit whose medians are in posteriors_invivo_withDIP.csv.
