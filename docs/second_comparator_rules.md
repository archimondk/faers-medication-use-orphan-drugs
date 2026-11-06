# Second specialty-pharmacy comparator set: selection rules and analysis plan

Version 1.0, 2026-11-06. This file is committed to the repository before any outcome data are produced for the drugs below. Its commit time is the record that the rules were fixed in advance. Any later change is made in a new commit with the reason stated.

Criterion C3 was worded after the FDA OOPD search in section 5. That search uses regulatory data only; no FAERS data of any kind had been produced for the candidate drugs when C3 was fixed.

## 1. Purpose

The existing specialty comparator set contains five drugs: teriflunomide, dimethyl fumarate, fingolimod, apremilast and tofacitinib. All five treat multiple sclerosis or immune-mediated disease. A second channel-matched set from a different therapeutic area tests whether the rare-versus-specialty estimates depend on that disease area.

## 2. Candidate universe

The candidates were defined by therapeutic class, not by individual drug:

- androgen-receptor pathway inhibitors for prostate cancer, and
- CDK4/6 inhibitors for hormone-receptor-positive breast cancer.

These classes were chosen for four reasons:

- They are long-term oral therapies for two common cancers.
- They are dispensed through specialty pharmacies.
- Each class has at least three branded agents from different manufacturers.
- They lie outside multiple sclerosis and immune-mediated disease.

Other oral oncology classes were considered and not used:

- HER2 tyrosine kinase inhibitors: neratinib and tucatinib have FDA orphan designations for breast cancer brain metastases [S9, S10].
- PARP inhibitors: these are partly indicated for ovarian cancer, and olaparib holds an orphan designation for ovarian cancer at least in the EU [S11].
- PI3K/AKT inhibitors: alpelisib is also marketed for an orphan indication (Vijoice), which makes name-based identification ambiguous. Capivasertib was approved in 2023.

## 3. Inclusion criteria

A drug in the candidate universe is included if it meets all of the following criteria. Each was assessed from regulatory and public sources, without FAERS outcome data.

| # | Criterion | Source |
|---|---|---|
| C1 | Listed on at least one of the two specialty drug lists already used for Table S1b (L1 Accredo/Express Scripts 2025, L2 CVS Caremark April 2026) | L1, L2 |
| C2 | Oral and self-administered, for long-term use | US prescribing information |
| C3 | The principal US indication is not orphan-designated. Designations for other indications are recorded but do not exclude the drug. These include unapproved indications, combination products and withdrawn designations. | FDA OOPD database [S0], manual search (section 5); review of CDK inhibitors [S8] |
| C4 | No REMS with restricted distribution | US prescribing information, FDA REMS list |
| C5 | No US generic marketed before 2025Q4, so the distribution channel did not shift to retail during the study period | Orange Book, ANDA approval letters, generic-availability listings |
| C6 | US approval no later than 31 December 2020, giving at least five years of reports in the study period | Drugs@FDA |
| C7 | Unambiguous identification by generic and brand name in DRUGNAME/PROD_AI | Count screen, file S3 |
| C8 | Manufacturer patient support programme documented | Manufacturer pages |

**Dosing regimen.** Drugs taken continuously (every day) form the main set. Drugs given on a 21-days-on, 7-days-off cycle are analysed only as a sensitivity addition, because the scheduled off-week may be reported as an omitted dose.

## 4. Assessment of candidates

| Drug (brand; manufacturer) | US approval | Regimen | C1 L1 / L2 | C3 orphan | C4 REMS | C5 generic | C8 programme | Assignment |
|---|---|---|---|---|---|---|---|---|
| Enzalutamide (Xtandi; Astellas, co-marketed by Pfizer) | Aug 2012 | 160 mg once daily, continuous | Yes / Yes | None (OOPD) | No | ANDAs approved 2021–2026, none marketed; Sandoz final approval 20 Apr 2026 [S1, S2] | Xtandi Support Solutions [S12] | Main set |
| Apalutamide (Erleada; Janssen) | Feb 2018 | 240 mg once daily, continuous | Yes / Yes | None (OOPD) | No | Zydus generic approved 17 Mar 2025, not commercially available [S3] | Janssen CarePath [S13] | Main set |
| Darolutamide (Nubeqa; Bayer) | Jul 2019 | 600 mg twice daily, continuous | Yes / Yes | None (OOPD) | No | None approved [S4] | Access Services by Bayer [S14] | Main set |
| Abemaciclib (Verzenio; Lilly) | Sep 2017 | Twice daily, continuous | Yes / No | None (OOPD) | No | None approved [S5] | Lilly Oncology Support Center [S15] | Main set |
| Palbociclib (Ibrance; Pfizer) | Feb 2015 | 125 mg once daily, 21 days on / 7 off | Yes / Yes | Chordoma, 2020, for the combination afatinib + palbociclib; not approved (OOPD) | No | None approved [S6] | Pfizer Oncology Together [S16] | Cyclic sensitivity set |
| Ribociclib (Kisqali; Novartis) | Mar 2017 | Once daily, 21 days on / 7 off | Yes / Yes | Neuroblastoma, 2022; not approved (OOPD) | No | None approved [S7] | Kisqali patient support [S17] | Cyclic sensitivity set |
| Abiraterone (Zytiga, Yonsa; Janssen, Sun) | Apr 2011 | Once daily, continuous | Yes / generic only | — | No | US generics launched 23 Nov 2018 [S18] | — | Excluded (C5) |

Notes:

- The Kisqali Femara Co-Pack is identified as ribociclib.
- Pfizer markets palbociclib and co-markets enzalutamide, so no manufacturer has more than one drug in the main set.
- The main set has four drugs from four manufacturers.

## 5. Checks completed before the analysis run

**FDA OOPD search**

- Search details:
  - Date: 2026-11-06.
  - Searched by: Xia Yuhang.
  - Environment: Windows 11, Chrome, Ego Lite 0.6.2.7.
  - Database: https://www.accessdata.fda.gov/scripts/opdlisting/oopd/index.cfm.
  - Settings: all designations, 1983-01-01 to 2026-11-06, condensed list.
  - Method: each drug was searched by generic name and by brand name.
- Record identifiers are the cfgridkey values of the OOPD detail pages, because the public pages do not show a designation number.

| Set | Drug | Designations found (cfgridkey) | Status | Assessment under C3 |
|---|---|---|---|---|
| Second, main | Enzalutamide | None | — | Meets |
| Second, main | Apalutamide | None | — | Meets |
| Second, main | Darolutamide | None | — | Meets |
| Second, main | Abemaciclib | None | — | Meets |
| Second, cyclic | Palbociclib | Treatment of chordoma, 2020-11-12, combination "afatanib and palbociclib" (781420) | Designated | Meets; combination product, unapproved indication |
| Second, cyclic | Ribociclib | Treatment of neuroblastoma, 2022-11-25 (912722) | Designated | Meets; unapproved indication |
| First | Teriflunomide | None | — | Meets |
| First | Dimethyl fumarate | Treatment of Friedreich's ataxia, 2013-09-11 (406713); MS in patients with flushing, fixed combination with aspirin 162 mg, 2018-09-27 (607917) | Designated | Meets; unapproved indication and combination product |
| First | Fingolimod | Treatment of chronic inflammatory demyelinating polyneuropathy, 2010-04-30 (305910) | Withdrawn or revoked | Meets |
| First | Apremilast | Treatment of Behçet's disease, 2013-01-17 (384812); pediatric ulcerative colitis, 2018-01-17 (622417) | Behçet's: designated and approved; pediatric UC: designated | Meets; principal indication psoriasis and psoriatic arthritis, but holds one approved orphan indication |
| First | Tofacitinib | None | — | Meets |
| Excluded (C5) | Abiraterone | None | — | Not applicable |

Consequences:

- A rule excluding any designation would have removed dimethyl fumarate, fingolimod and apremilast from the first set and palbociclib and ribociclib from the second. C3 is therefore defined on the principal indication. The stricter definition is kept as a sensitivity analysis (section 7).
- Apremilast is the only comparator in either set with an approved orphan indication. Its influence is already examined by the leave-one-out analysis of the first set.
- The description of the first set as "non-orphan" in the manuscript is replaced by a statement that its drugs are used mainly for common diseases and that three hold orphan designations for other indications.

**Generic entry in the first set**

Recorded 2026-11-06 from public sources, before any outcome data were produced for the configurations below. The sensitivity analysis drops each drug's reports from the FAERS quarterly file of the quarter with the earliest evidence of US generic marketing onward (configuration `configs/c5_specialty1_truncated.json`).

| Drug | First US generic marketing | Reports dropped from | Source |
|---|---|---|---|
| Dimethyl fumarate | Mylan launch, 19 August 2020 | 2020Q3 | S19 |
| Fingolimod | Generics notified 29 September 2022, confirmed launched by Novartis on 25 October 2022 | 2022Q3 | S20, S21 |
| Teriflunomide | Sandoz and Accord launches confirmed 14 March 2023 | 2023Q1 | S22 |
| Apremilast | None; generics barred until 2028 by patent ruling | Not truncated | S23 |
| Tofacitinib | Generic tablets first marketed in May 2026, after the study period | Not truncated | S24 |

Dimethyl fumarate, fingolimod and teriflunomide do not meet C5 as written. This is reported as a limitation of the first set (see section 7).

## 6. Count screen (no outcome data)

`scripts/screen_comparator_counts.py` reads only DEMO and DRUG files. For each candidate, and for the existing groups for reference, it reports:

- de-duplicated primary-suspect reports and their distribution by year;
- report type, reporter and country mix;
- the drug-name strings matched.

Decisions allowed after the screen, and no others:

- A drug with fewer than 100 primary-suspect reports in 2013Q1–2025Q4 is kept at report level but left out of drug-level analyses, the same threshold used for the rare-disease drugs.
- Name patterns are corrected if S3 shows matches to unrelated products. The corrected patterns are committed before the analysis run.
- No drug is excluded for its report mix or its volume.

### Screen results (run 2026-11-06, 2013Q1–2025Q4)

**Check against the main analysis.** The screen reproduced the report totals of Table 2 exactly, together with their consumer and periodic shares:

- rare-disease group: 15,508 reports (62.0% consumer, 60.5% periodic);
- first specialty set: 488,215 reports (49.7% consumer, 68.0% periodic);
- common comparators: 174,358 reports (41.5% consumer, 23.3% periodic).

**Volume.** Every candidate had at least 100 primary-suspect reports:

| Drug | Reports | Share of the main second set |
|---|---|---|
| Enzalutamide | 47,286 | 63.2% |
| Apalutamide | 9,142 | 12.2% |
| Darolutamide | 3,059 | 4.1% |
| Abemaciclib | 15,313 | 20.5% |
| Palbociclib | 80,180 | cyclic set |
| Ribociclib | 20,275 | cyclic set |

**Name patterns.** The matched strings were the US brand and generic names, development codes (JNJ-56021927, ARN-509, LEE011, PD-0332991) and non-US brands (Verzenios, Ramiven, Palbace, Reampla, Etulsi, Kryxana, Valamor, Risarg). Every match carried the expected PROD_AI. Blinded trial entries were few: "BLINDED ENZALUTAMIDE" (74 reports), plus small numbers for palbociclib, ribociclib and abemaciclib. They are kept, as they are for every group in the main analysis. No report matched more than one screened drug. The patterns are unchanged.

**Decisions.** No decision under this section was triggered. The rules and patterns are unchanged.

**Report mix recorded before the outcome run.** These observations are recorded here so that they cannot be chosen after the results are known:

| Group | Reports | Expedited % | Periodic % | Consumer % | Documented HCP % | HP/RN % | US % |
|---|---|---|---|---|---|---|---|
| Rare-disease | 15,508 | 37.4 | 60.5 | 62.0 | 25.8 | 8.3 | 81.7 |
| First specialty set | 488,215 | 26.4 | 68.0 | 49.7 | 34.4 | 14.2 | 84.9 |
| Second set, main | 74,800 | 49.1 | 42.2 | 63.5 | 28.1 | 6.1 | 73.3 |
| Second set + cyclic | 175,255 | 45.4 | 43.7 | 52.7 | 31.9 | 12.9 | 71.3 |
| Common | 174,358 | 71.3 | 23.3 | 41.5 | 37.9 | 10.3 | 42.8 |

- The main second set has fewer periodic reports than the rare-disease group or the first set, and more expedited reports. A difference between the two specialty-set estimates may therefore partly reflect report type. The adjusted model and the subset excluding periodic reports, both already in the plan, are reported next to the crude estimate.
- Enzalutamide provides 63% of the main second set. Of its reports, 74% were received in 2013–2018, followed by a sharp fall. The drug-level analysis, which is the primary comparison, and the leave-one-out analysis within the second set address its weight.

## 7. Analysis plan

All analyses use the main-analysis definitions unchanged:

- study period 2013Q1–2025Q4;
- de-duplication and primary-suspect attribution as in the main analysis;
- the reviewed PT dictionary and SMQ sensitivity definitions;
- drug-level REML with Hartung–Knapp;
- CR1 cluster-robust logistic regression with the same covariates, and t critical values with the number of drug clusters minus 1.

The following analyses are run with the main second set (four drugs) as the reference:

1. Report-level ROR for nonadherence, composite, dosing or administration error, and access or supply problems: rare-disease versus second set.
2. Drug-level meta-regression: rare-disease versus second set.
3. Adjusted report-level model with the second set as the reference group.
4. The Table 5 report subsets, using the same scenarios as for the first set.
5. The Nefecon worked example over the same three periods, with the second set's reference restricted to the same quarters.
6. Leave-one-out analysis within the second set.

Additional comparisons:

- Second set versus first set, at report and drug level. This compares the two channel-matched comparators directly.
- The pooled specialty reference (first set plus main second set, nine drugs), at report and drug level.
- The cyclic sensitivity set (main second set plus palbociclib and ribociclib).
- The first set truncated at the first US generic entry of each drug, as a sensitivity analysis for C5.
- A strict non-orphan reference made up of the specialty drugs without any OOPD designation: teriflunomide, tofacitinib and the four drugs of the main second set. This is a sensitivity analysis for C3.

The ambiguous-interruption category stays separate from nonadherence. This matters for oncology drugs, whose doses are often interrupted for toxicity.

## 8. Prespecified interpretation

The primary quantity is the drug-level nonadherence OR for rare-disease drugs versus the second set (OR₂), compared with the same OR versus the first set (OR₁ = 1.21).

- If OR₂/OR₁ lies between 0.67 and 1.5, the attenuation seen with the first set is reported as reproduced with channel-matched drugs from a different therapeutic area.
- If OR₂/OR₁ lies outside that range, the result is reported as showing that the estimate depends on which channel-matched drugs are chosen. This supports the main conclusion that no single comparator gives a stable interpretation.
- In neither case is the result described as establishing the source of the reports, such as patient support programmes.

## Sources

- S1. FDA. ANDA 216068 approval letter, enzalutamide tablets (Sandoz), 20 April 2026. https://www.accessdata.fda.gov/drugsatfda_docs/appletter/2026/216068Orig1s000Panoramaltr.pdf
- S2. Drugs.com. Generic Xtandi availability. https://Drugs.com/availability/generic-xtandi.html; OncLive, 28 June 2026, tentative approval of Lupin enzalutamide tablets. https://www.onclive.com/view/generic-enzalutamide-tablets-granted-tentative-fda-approval-for-prostate-cancer
- S3. Drugs.com. Generic Erleada availability. https://drugs.com/availability/generic-erleada.html
- S4. Drugs.com. Generic Nubeqa availability. https://drugs.com/availability/generic-nubeqa.html
- S5. Drugs.com. Generic Verzenio availability. https://Drugs.com/availability/generic-verzenio.html
- S6. Drugs.com. Generic Ibrance availability. https://Drugs.com/availability/generic-ibrance.html
- S7. Drugs.com. Generic Kisqali availability. https://drugs.com/availability/generic-kisqali.html
- S0. FDA Office of Orphan Products Development. Orphan drug designations and approvals database, searched 2026-11-06. https://www.accessdata.fda.gov/scripts/opdlisting/oopd/index.cfm
- S8. CDK inhibitors and FDA: approved and orphan. Cancers. 2024;16(8):1555. https://www.mdpi.com/2072-6694/16/8/1555
- S9. Healio. FDA grants orphan drug designation to tucatinib for brain metastases. 2017. https://www.healio.com/news/hematology-oncology/20170613/fda-grants-orphan-drug-designation-to-tucatinib-for-brain-metastases
- S10. OncLive. FDA grants neratinib orphan drug status for breast cancer with brain metastases. https://www.onclive.com/view/fda-grants-neratinib-orphan-drug-status-for-breast-cancer-with-brain-metastases
- S11. EMA. EU/3/07/501 orphan designation for treatment of ovarian cancer (olaparib). https://www.ema.europa.eu/en/medicines/human/orphan-designations/eu307501
- S12. Astellas. Xtandi Support Solutions. https://www.xtandi.com/pdf/xtandi-support-solutions-handout.pdf
- S13. Janssen. Starting and staying on Erleada resource guide. https://erleadahcp.com/documents/starting_and_staying_on_erleada_resource_guide.pdf
- S14. Bayer. Savings and support with Access Services by Bayer (Nubeqa). https://nubeqa-us.com/savings_and_support_with_access_services_by_bayer.pdf
- S15. Lilly. Verzenio specialty pharmacy and support. https://verzenio.lilly.com/savings-support/specialty-pharmacy
- S16. Pfizer. Ibrance personalized patient support (Pfizer Oncology Together). https://ibrance.pfizerpro.com/support-order/personalized-patient-support
- S17. Novartis. Kisqali patient support. https://support.kisqali.com
- S18. OptumRx. New generics: Zytiga, 27 November 2018. https://prod-prorx.optum.com/content/dam/optum3/professional-optumrx/news/rxnews/new-generics/newgenerics_zytiga_2018-1127.pdf
- S19. Drug Topics. Mylan launches generic version of multiple sclerosis treatment. 20 August 2020. https://www.drugtopics.com/view/mylan-launches-generic-version-of-multiple-sclerosis-treatment
- S20. OptumRx. New generics: Gilenya, 29 September 2022. https://professionals.optumrx.com/content/dam/optum3/professional-optumrx/news/rxnews/new-generics/newgenerics_gilenya_2022-0929.pdf
- S21. Generics Bulletin. ANDA sponsors silent as Novartis confirms US Gilenya generic entries. 28 October 2022. https://insights.citeline.com/GB152254/ANDA-Sponsors-Silent-As-Novartis-Confirms-US-Gilenya-Generic-Entries
- S22. Generics Bulletin. Sandoz, Accord confirm first US Aubagio generics as settlement day arrives. 14 March 2023. https://insights.citeline.com/GB152587/Sandoz-Accord-Confirm-First-US-Aubagio-Generics-As-Settlement-Day-Arrives
- S23. Amgen. Amgen wins patent appeal on Otezla (apremilast). April 2023. https://www-ext.amgen.com/newsroom/press-releases/2023/04/amgen-wins-patent-appeal-on-otezla-apremilast
- S24. Drugs.com. Generic Xeljanz availability. https://drugs.com/availability/generic-xeljanz.html
- L1. Accredo Health Group (Express Scripts). Specialty drug list. 2025. https://hr.duke.edu/sites/default/files/Specialty-Drug-List-2026.pdf
- L2. CVS Caremark. Advanced Control Specialty Formulary, April 2026. https://dhr.delaware.gov/benefits/cvs/documents/preferred-specialty-drug-list-040126.pdf
