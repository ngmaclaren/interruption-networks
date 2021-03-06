log using "Z:\Neil Interruptions\analysis-2SLS-extended-v5.smcl", replace

ssc install ranktest
ssc install ivreg2
ssc install ivreg2h
ssc install estout
ssc install coefplot
ssc install outreg2
ssc install mkcorr
ssc install rangestat

*import delimited data\speakingTime-data.csv, numericcols(3) clear
import delimited "Z:\Neil Interruptions\speakingTime-data-extended.csv", ///
numericcols(3) clear

gen d_male = 0
replace d_male = 1 if gender=="male"
gen d_english = 0
replace d_english = 1 if english_second_language=="yes"
gen d_simulation = 0
replace d_simulation = 1 if simulation=="cs"
gen d_institution = 0
replace d_institution = 1 if institution=="S2"
gen d_operator = 0
replace d_operator = 1 if participant_is_operator=="yes"

gen quiz = game_knowledge_quiz
*gen conscientiousness = conscientiousness
gen groupid = group_id
gen speaking = total_speaking_time/1000

*** Lewbel Estimation, w/ Interruptions Replacing Speaking Time ***

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss_count = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(iss_count) gen(hx)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss_count = intelligence d_operator conscientiousness-extraversion hx_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(orig_pr = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(orig_pr) gen(hy)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(orig_pr = intelligence d_operator conscientiousness-extraversion hy_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(new_pr = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(new_pr) gen(hz)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(new_pr = intelligence d_operator conscientiousness-extraversion hz_*), ///
cluster(groupid)
estat firststage
estat endogenous

*** Lewbel Estimation w/ Two Endog Variables ***
ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss_count speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(iss_count speaking) gen(ha)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(iss_count speaking = intelligence d_operator conscientiousness-extraversion ha_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(orig_pr speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(orig_pr speaking) gen(hb)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(orig_pr speaking = intelligence d_operator conscientiousness-extraversion hb_*), ///
cluster(groupid)
estat firststage
estat endogenous

ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(new_pr speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(new_pr speaking) gen(hc)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(new_pr speaking = intelligence d_operator conscientiousness-extraversion hc_*), ///
cluster(groupid)
estat firststage
estat endogenous


log close
translate "Z:\Neil Interruptions\analysis-2SLS-extended-v5.smcl" "Z:\Neil Interruptions\analysis-2SLS-extended-v5.log", replace linesize(100) translator(smcl2log)


*** Linear Analysis ***
* demonstrate weak instruments
ivreg2 planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(speaking)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid)
estat firststage
estat endogenous

* the model with Lewbel instruments
ivreg2h planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator conscientiousness-extraversion), ///
cluster(groupid) first endog(speaking) gen(h)

separate planning_phase_vote_total, by(d_male)

quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator conscientiousness-extraversion h_*), ///
cluster(groupid)
estat firststage
estat endogenous
 
quietly margins, over(d_male) at(speaking=(0(10)280))
marginsplot, ///
  addplot(scatter planning_phase_vote_total0 speaking, ///
     msize(small) msymbol(circle) ///
     || scatter planning_phase_vote_total1 speaking, ///
     msize(small) msymbol(diamond)) ///
  title("") ///
  xtitle("Speaking Time (seconds)") ytitle("Leader Emergence (votes)") ///
  graphregion(fcolor(white)lcolor(white)) scheme(s2mono) ///
  plotdimension( , labels("Females, margins estimate" "Males, margins estimate"))

** exploratory: same regression with operator as an included instrument **
ivreg2h planning_phase_vote_total ///
d_operator d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence conscientiousness-extraversion), ///
cluster(groupid) first endog(speaking)
* d_operator highly significant in the first stage and not significant in the
* second stage

** same for intelligence
ivreg2h planning_phase_vote_total ///
intelligence d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = d_operator consc-extraversion), ///
cluster(groupid) first endog(speaking)
* same for intelligence

** same for openness
ivreg2h planning_phase_vote_total ///
openness d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator consc-neuroticism extraversion), ///
cluster(groupid) first endog(speaking)
* same for openness

*** Post Hoc Analysis ***
sem (planning_phase_vote_total <- speaking ///
d_male age quiz d_english group_size d_simulation d_institution) ///
(speaking <- intelligence d_operator conscientiousness-extraversion h_* ///
d_male age quiz d_english group_size d_simulation d_institution), ///
vce(cluster groupid) cov(e.planning_phase_vote_total*e.speaking)

test _b[/cov(e.planning_phase_vote_total,e.speaking)] = 0

test ///
(_b[speaking: intelligence ] = 0) ///
(_b[speaking: d_operator ] = 0) ///
(_b[speaking: conscientiousness ] = 0) ///
(_b[speaking: agreeableness ] = 0) ///
(_b[speaking: neuroticism ] = 0) ///
(_b[speaking: openness ] = 0) ///
(_b[speaking: extraversion ] = 0) ///
(_b[speaking: h_speaking_d_male_g ] = 0) ///
(_b[speaking: h_speaking_age_g ] = 0) ///
(_b[speaking: h_speaking_quiz_g ] = 0) ///
(_b[speaking: h_speaking_d_english_g ] = 0) ///
(_b[speaking: h_speaking_group_size_g ] = 0) ///
(_b[speaking: h_speaking_d_simulation_g ] = 0) ///
(_b[speaking: h_speaking_d_institution_g ] = 0)

dis "F-test equivalent = "  r(chi2)/r(df)

estat gof, stats(all)

estat stdize : testnl ///
_b[planning_phase_vote_total:speaking]*_b[speaking:d_operator] = ///
_b[planning_phase_vote_total:speaking]*_b[speaking:openness] = ///
_b[planning_phase_vote_total:speaking]*_b[speaking:d_male] = ///
_b[planning_phase_vote_total:speaking]*_b[speaking:intelligence]

estat stdize: nlcom _b[planning_phase_vote_total:d_male] +  ///
_b[planning_phase_vote_total:speaking]* _b[speaking:d_male]

estat stdize: nlcom _b[planning_phase_vote_total:speaking]* _b[speaking:d_male]

** need to include John's "social skills" and "organizational skills" tests
* indirect effects to compare are:
* operator status
estat stdize: ///
nlcom _b[planning_phase_vote_total:speaking]*_b[speaking:d_operator]
* intelligence and openness
estat stdize: ///
nlcom _b[planning_phase_vote_total:speaking]*_b[speaking:intelligence]
estat stdize: ///
nlcom _b[planning_phase_vote_total:speaking]*_b[speaking:openness]
* extraversion, agreeableness, neuroticism
estat stdize: ///
nlcom _b[planning_phase_vote_total:speaking]*_b[speaking:extraversion]
estat stdize: ///
nlcom _b[planning_phase_vote_total:speaking]*_b[speaking:agreeableness]
estat stdize: ///
nlcom _b[planning_phase_vote_total:speaking]*_b[speaking:neuroticism]
* conscientiousness
estat stdize: ///
nlcom _b[planning_phase_vote_total:speaking]*_b[speaking:conscientiousness]

estat stdize: ///
testnl _b[planning_phase_vote_total:speaking]*_b[speaking:d_operator] = /// 
(_b[planning_phase_vote_total:speaking]*_b[speaking:intelligence] ///
+ _b[planning_phase_vote_total:speaking]*_b[speaking:openness])

*** For Appendix ***
** reduced model with cluster-robust DWH
ivreg2 planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator openness), ///
cluster(groupid) first endog(speaking)
quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator openness), ///
cluster(groupid) first
* outreg2 using appendix.doc, replace ctitle(Reduced Model)
estat firststage
estat endogenous

** reduced model (no openness) with cluster-robust DWH
ivreg2 planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator), ///
cluster(groupid) first endog(speaking)
quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator), ///
cluster(groupid) first
* outreg2 using appendix.doc, replace ctitle(Reduced Model)
estat firststage
estat endogenous

** original model, but with remaining personality vars as included inst
ivreg2 planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
consc-neuroticism extraversion ///
(speaking = intelligence d_operator openness), ///
cluster(groupid) first endog(speaking)
quietly ivregress 2sls planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
consc-neuroticism extraversion ///
(speaking = intelligence d_operator openness), ///
cluster(groupid)
estat firststage
estat endogenous

** Poisson model with cluster-robust DWH
ivpoisson cfunction planning_phase_vote_total ///
d_male age quiz d_english group_size d_simulation d_institution ///
(speaking = intelligence d_operator consc-extraversion h_*), ///
vce(cluster groupid)
* outreg2 using appendix.doc, append ctitle(Poisson Model)
test _b[/c_speaking] = 0

*** Summary Statistics and Plots ***
summarize planning_phase_vote_total speaking ///
d_male age quiz d_english group_size d_simulation d_institution ///
intelligence d_operator conscientiousness-extraversion h_*

correlate planning_phase_vote_total speaking ///
d_male age quiz d_english group_size d_simulation d_institution ///
intelligence d_operator conscientiousness-extraversion h_*

*mkcorr planning_phase_vote_total speaking ///
*d_male age quiz d_english group_size d_simulation d_institution ///
*intelligence d_operator consc-extraversion h_*, ///
*log(summaryStatistics) replace means num cdec(2) mdec(2) casewise

loneway speaking groupid
loneway planning_phase_vote_total groupid

kdensity planning_phase_vote_total, ///
graphregion(fcolor(white)lcolor(white)) scheme(s2mono) normal
* graph export kdensityVotes.png

kdensity speaking, ///
graphregion(fcolor(white)lcolor(white)) scheme(s2mono) normal
* graph export kdensitySpeaking.png

scatter planning_phase_vote_total speaking, ///
graphregion(fcolor(white)lcolor(white)) scheme(s2mono) ///
|| lowess planning_phase_vote_total speaking ///
|| lfit planning_phase_vote_total speaking
* graph export scatter.png

centile speaking, centile(61.67)
pcorr planning_phase_vote_total speaking d_male

