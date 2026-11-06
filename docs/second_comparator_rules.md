# Second specialty-pharmacy comparator set: selection rules and analysis plan

Version 1.0, 2026-11-06. This file is committed to the repository before any outcome data are produced for the drugs below. Its commit time is the record that the rules were fixed in advance. Any later change is made in a new commit with the reason stated.

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
| C3 | No FDA orphan designation for any indication | FDA OOPD database (manual search, see section 5) |
| C4 | No REMS with restricted distribution | US prescribing information, FDA REMS list |
| C5 | No US generic marketed before 2025Q4, so the distribution channel did not shift to retail during the study period | Orange Book, ANDA approval letters, generic-availability listings |
| C6 | US approval no later than 31 December 2020, giving at least five years of reports in the study period | Drugs@FDA |
| C7 | Unambiguous identification by generic and brand name in DRUGNAME/PROD_AI | Count screen, file S3 |
| C8 | Manufacturer patient support programme documented | Manufacturer pages |

**Dosing regimen.** Drugs taken continuously (every day) form the main set. Drugs given on a 21-days-on, 7-days-off cycle are analysed only as a sensitivity addition, because the scheduled off-week may be reported as an omitted dose.

## 4. Assessment of candidates

| Drug (brand; manufacturer) | US approval | Regimen | C1 L1 / L2 | C3 orphan | C4 REMS | C5 generic | C8 programme | Assignment |
|---|---|---|---|---|---|---|---|---|
| Enzalutamide (Xtandi; Astellas, co-marketed by Pfizer) | Aug 2012 | 160 mg once daily, continuous | Yes / Yes | None found; OOPD to confirm | No | ANDAs approved 2021–2026, none marketed; Sandoz final approval 20 Apr 2026 [S1, S2] | Xtandi Support Solutions [S12] | Main set |
| Apalutamide (Erleada; Janssen) | Feb 2018 | 240 mg once daily, continuous | Yes / Yes | None found; OOPD to confirm | No | Zydus generic approved 17 Mar 2025, not commercially available [S3] | Janssen CarePath [S13] | Main set |
| Darolutamide (Nubeqa; Bayer) | Jul 2019 | 600 mg twice daily, continuous | Yes / Yes | None found; OOPD to confirm | No | None approved [S4] | Access Services by Bayer [S14] | Main set |
| Abemaciclib (Verzenio; Lilly) | Sep 2017 | Twice daily, continuous | Yes / No | None found [S8]; OOPD to confirm | No | None approved [S5] | Lilly Oncology Support Center [S15] | Main set |
| Palbociclib (Ibrance; Pfizer) | Feb 2015 | 125 mg once daily, 21 days on / 7 off | Yes / Yes | None found [S8]; OOPD to confirm | No | None approved [S6] | Pfizer Oncology Together [S16] | Cyclic sensitivity set |
| Ribociclib (Kisqali; Novartis) | Mar 2017 | Once daily, 21 days on / 7 off | Yes / Yes | None found [S8]; OOPD to confirm | No | None approved [S7] | Kisqali patient support [S17] | Cyclic sensitivity set |
| Abiraterone (Zytiga, Yonsa; Janssen, Sun) | Apr 2011 | Once daily, continuous | Yes / generic only | — | No | US generics launched 23 Nov 2018 [S18] | — | Excluded (C5) |

Notes:

- The Kisqali Femara Co-Pack is identified as ribociclib.
- Pfizer markets palbociclib and co-markets enzalutamide, so no manufacturer has more than one drug in the main set.
- The main set has four drugs from four manufacturers.

## 5. Checks still to be completed before the analysis run

1. **OOPD.** Search https://www.accessdata.fda.gov/scripts/opdlisting/oopd/ by generic name for each of the six candidates. Record the search date and any designation found. A designation for any indication excludes the drug under C3. The secondary sources above found none, but OOPD has not been searched directly.
2. **Consistency with the first set.** Run the same OOPD search for the five existing specialty comparators. This applies C3 identically to both sets.
3. **Generic entry in the first set.** Record when US generics were first marketed for the existing specialty comparators. Dimethyl fumarate (about 2020), fingolimod (about 2022) and teriflunomide (about 2023) are to be verified. These drugs do not meet C5 as written, and this is reported as a limitation of the first set (see section 7).

## 6. Count screen (no outcome data)

`scripts/screen_comparator_counts.py` reads only DEMO and DRUG files. For each candidate, and for the existing groups for reference, it reports:

- de-duplicated primary-suspect reports and their distribution by year;
- report type, reporter and country mix;
- the drug-name strings matched.

Decisions allowed after the screen, and no others:

- A drug with fewer than 100 primary-suspect reports in 2013Q1–2025Q4 is kept at report level but left out of drug-level analyses, the same threshold used for the rare-disease drugs.
- Name patterns are corrected if S3 shows matches to unrelated products. The corrected patterns are committed before the analysis run.
- No drug is excluded for its report mix or its volume.

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
- L1. Accredo Health Group (Express Scripts). Specialty drug list. 2025. https://hr.duke.edu/sites/default/files/Specialty-Drug-List-2026.pdf
- L2. CVS Caremark. Advanced Control Specialty Formulary, April 2026. https://dhr.delaware.gov/benefits/cvs/documents/preferred-specialty-drug-list-040126.pdf
