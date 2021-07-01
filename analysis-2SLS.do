log using "analysis-2SLS.smcl", replace

ssc install ranktest
ssc install ivreg2
ssc install ivreg2h
ssc install estout
ssc install coefplot
ssc install outreg2
ssc install mkcorr
ssc install rangestat

import delimited "speakingTime-data-extended.csv", numericcols(3) clear

gen quiz = game_knowledge_quiz
gen groupid = group_id
gen speaking = total_speaking_time

*** Lewbel Estimation, w/ Interruptions Replacing Speaking Time ***

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(iss) gen(hx)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss = intelligence d_operator conscientiousness-extraversion hx_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(nss = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(nss) gen(hy)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(nss = intelligence d_operator conscientiousness-extraversion hy_*), ///
cluster(groupid)
estat firststage
estat endogenous

* i_pr, n_pr, b_pr
ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(i_pr = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(i_pr) gen(hi)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(i_pr = intelligence d_operator conscientiousness-extraversion hi_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(n_pr = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(n_pr) gen(hn)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(n_pr = intelligence d_operator conscientiousness-extraversion hn_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_pr = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(b_pr) gen(hb)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_pr = intelligence d_operator conscientiousness-extraversion hb_*), ///
cluster(groupid)
estat firststage
estat endogenous

*** Lewbel Estimation w/ Two Endog Variables ***
ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(iss speaking) gen(ha)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss speaking = intelligence d_operator conscientiousness-extraversion ha_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(nss speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(nss speaking) gen(hb)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(nss speaking = intelligence d_operator conscientiousness-extraversion hb_*), ///
cluster(groupid)
estat firststage
estat endogenous

* i_pr, n_pr, b_pr
ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(i_pr speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(i_pr speaking) gen(his)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(i_pr speaking = intelligence d_operator conscientiousness-extraversion his_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(n_pr speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(n_pr speaking) gen(hns)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(n_pr speaking = intelligence d_operator conscientiousness-extraversion hns_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_pr speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(b_pr speaking) gen(hbs)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_pr speaking = intelligence d_operator conscientiousness-extraversion hbs_*), ///
cluster(groupid)
estat firststage
estat endogenous

* demonstrating that using weighted in-degree doesn't work
ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_wid speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(b_wid speaking) gen(hwids)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_wid speaking = intelligence d_operator conscientiousness-extraversion hwids_*), ///
cluster(groupid)
estat firststage
estat endogenous

* demonstrating that the p-value for the PageRank coefficient goes down if 
* those values are considered zero for unconnected nodes (e.g., if NA -> 0)
replace b_pr = 0 if b_wid == 0 & b_wod == 0

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_pr speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(b_pr speaking) gen(hb2s)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(b_pr speaking = intelligence d_operator conscientiousness-extraversion hb2s_*), ///
cluster(groupid)
estat firststage
estat endogenous

log close
translate "analysis-2SLS.smcl" "analysis-2SLS.log", replace linesize(100) translator(smcl2log)