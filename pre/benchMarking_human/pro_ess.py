import pandas as pd

ess_lnc = pd.read_csv("../../data/benchMarking/esslnc.csv")


human_ess_lnc = ess_lnc[ess_lnc['Organism'] == 'Human']
mouse_ess_lnc = ess_lnc[ess_lnc['Organism'] == 'Mouse']

human_ess_lnc = human_ess_lnc[(human_ess_lnc['cancer_related'] > 0)]
mouse_ess_lnc = mouse_ess_lnc[(mouse_ess_lnc['cancer_related'] > 0)]

human_ess_lnc = human_ess_lnc[['Noncode_id','Lncbook_id','lib_id','gene_name','ensembl_id']]
mouse_ess_lnc = mouse_ess_lnc[['Noncode_id','Lncbook_id','lib_id','gene_name','ensembl_id']]

human_ess_lnc.to_csv("../../data/benchMarking/human/esslnc_lit.csv", index=False)
mouse_ess_lnc.to_csv("../../data/benchMarking/mouse/esslnc_lit.csv", index=False)