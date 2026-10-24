# Stage2 registration: V6 Stage2a

## Background (five lines)

Stage1 used 136 complete DEV participants (544 records) and frozen original PRRMrel B-cal readouts.
The primary deployed readout uses training z statistics and a ±3 clip; clipping affects short subrecords more strongly than full records.
The DEV decision, made before new CONF computations, assigns unclipped indices to variance/noise/agreement and clipped indices to SDS accuracy.
Stage1 showed across-record independent-error tension, while stored-position, adjacency and recorded-quality associations were not detected by the specified criteria.
Stage2a supplies DEV diagnostics and fixed numeric CONF predictions; Stage2b is locked until a manual registration receipt is supplied.

## Data, split and exposure

{"dataset_people": 276, "records": 1092, "DEV_people": 141, "DEV_complete": 136, "DEV_complete_records": 544, "CONF_people": 135, "CONF_complete": 131, "analysed_record": "130 seconds, samples1280:17920 at128Hz", "split_sha256": "69ece5f8828f7296374cc82125e0294ed8eef8cd5ff772f477cd2b0b8cb9606d", "position_definition": "stored session number (recording order unverified)"}

No new CONF features, predictions, thinning, quality variables, accuracy, intervals or statistics were computed in V5 Stage1 or V6 Stage2a. Earlier full-sample analyses already included CONF: four-record ICC, single-record r, four-record r and participant-mean thinning. Frozen fold readouts used CONF labels in training. CLEAN means the registered quantity was not previously computed full-sample; it does not mean a pristine independent holdout.

## Fixed computations

- **timestamp_utc**: 2026-10-24T11:47:21Z
- **stage**: V6 Stage2a DEV only
- **stage1_manifest_sha256**: 52d5927f9e4c0a87db8423bbe401f39f806592800bc8594cdf2d47ad627199d9
- **DEV_complete_n**: 136
- **CONF_complete_n**: 131
- **bootstrap_DEV**: exact Stage1 2000x136 matrix, seed20261024
- **bootstrap_CONF**: {'seed': 20261025, 'draws': 2000}
- **thinning_seed_base**: 20261101
- **thinning_seed_index**: zero-based index within sorted complete cohort (participant,session), reset separately for each cohort
- **thinning_rng**: numpy default_rng PCG64, channel-major then original cached frequency-major; one uniform per event; retain u<q; q=.8 nested in q=.9; same thinned trains for both readouts
- **full_accuracy**: rstar: subtract the mean of each stored session position, then correlate pooled predictions with repeated SDS; recenter within each bootstrap
- **subrecord_accuracy**: raw pooled H1/H2, Q1-Q4, C1/C2 or H13_1/H13_2 against repeated SDS; no additional centring; equal rows per participant
- **quarter_ICC**: records as rows and quarters as columns; two-way consistency ICC; bootstrap people with all their records/quarters
- **variance_primary**: unclipped; clipped sensitivity
- **SDS_accuracy_primary**: clipped; unclipped sensitivity
- **thinning_primary**: unclipped; clipped sensitivity
- **prediction_SE**: sample SD (ddof1) of finite 2000 bootstrap values
- **prediction_interval**: DEV estimate +/- z * SE * sqrt(1+136/131); z=1.645/1.96; no clipping to admissible parameter range
- **numeric_pass**: CONF point estimate inside closed 95% prediction interval
- **categorical_pass**: registered boolean conditions or exact DEV label match; no numeric prediction interval for a category
- **overall**: LOOCV only; 46fold reported separately, never selects the main verdict
- **missing**: no imputation; exclude only observations needed by a metric, report counts; undefined estimates cannot pass
- **no_CONF_computation**: True
- **no_Section_A_estimates_when_written**: True
- **position**: stored session number (recording order unverified)

Use the same frozen A2 masked computation and saved fold readouts. No fitting or SDS tuning occurs. All bootstrap resamples keep each sampled participant’s four records and all subrecords. The half difference variance uses ddof=1. Quarter consistency ICC has records as rows and quarters as columns. Quality measures retain the Stage1 native scale (observed DEV CQ0–4), interpolation fraction and counter-gap convention.

## Statistic definitions

For full predictions Y_ij, let m_j be the mean at stored session number j, person means be Ybar_i, and grand mean be Ybar. All means and ranks are recomputed within each participant bootstrap resample.

- rstar = corr(Y_ij - m_j, SDS_i), pooled over records; r4 = corr(mean_j Y_ij, SDS_i).
- r_H and r_Q are Pearson correlations of pooled, uncentred half/quarter predictions against repeated SDS. Each participant contributes8 half rows or16 quarter rows. Contiguous and13s halves follow the same pooling rule.
- r_HH = corr(H1,H2) across records; R_SB = 2*r_HH/(1+r_HH). sigma2_half = sample Var(H1-H2)/2. sigma2_quarter = mean across records of the sample variance of Q1-Q4 (ddof1).
- gamma = log(sigma2_quarter/sigma2_half)/log(mean N_eff_half / mean N_eff_quarter). N_eff is the frozen A2 sum of within-block window lengths/(tau+1) over source events and selected coordinates.
- MSR = 4*sample Var(person means); MSW = mean of within-person sample variances; MSE = sum((Y_ij-Ybar_i-m_j+Ybar)^2)/(3*(n-1)). ICC(1)=(MSR-MSW)/(MSR+3*MSW); ICC(C,1)=(MSR-MSE)/(MSR+3*MSE).
- sigma2_s_counting = sigma2_half * mean N_eff_half / mean N_eff_full; S_counting = sigma2_s_counting/MSE; Delta_state = R_SB - ICC(C,1). These raw values are not truncated.
- r_QQ uses the same two-way consistency ICC formula with records as rows and the4 quarters as columns. D_half=r_H^2-r_HH; D_quarter=r_Q^2-r_QQ.
- D1=rstar^2-ICC(C,1); D4=r4^2-4*ICC(C,1)/(1+3*ICC(C,1)).
- LMM compares y~C(stored position)+(1|participant) with y~1+(1|participant), by Gaussian ML, nonnegative random-intercept variance, likelihood-ratio test df3. Balanced orthogonal within/between ML is the Stage1 implementation.
- lag1-minus-lag3 is mean(corr(Y1,Y2),corr(Y2,Y3),corr(Y3,Y4))-corr(Y1,Y4), without person-centring.
- P5 comprises Spearman correlations of native mean contact quality, interpolation fraction and counter-gap count against squared two-way full-record residuals and squared H1-H2 differences.
- P8 uses position-centred rstar(q) and person-mean r4(q), divided by their corresponding unthinned estimates. The same thinning realization is shared by both readouts and both index variants.

## Prediction intervals and pass rules

SE is the sample SD (ddof=1) of the2,000 DEV bootstrap estimates. PI = theta_DEV ± z·SE·sqrt(1+136/131), z=1.645/1.96. Numeric predictions pass when the CONF point estimate is inside the closed95% PI. These are prediction intervals, not the DEV percentile confidence intervals. Categorical rows have no numeric PI; their exact pass rule is stated. All intervals are untruncated.

| ID | Definition | Readout / index | Tag | DEV estimate | Bootstrap SE | 90% prediction interval | 95% prediction interval | Pass rule |
|---|---|---|---|---:|---:|---|---|---|
| P1a | r_HH | LOOCV / unclipped | CLEAN | 0.012657992 | 0.041307951 | [-0.084352685, 0.10966867] | [-0.1029292, 0.12824518] | CONF estimate lies inside the closed 95% prediction interval |
| P1b | sigma2_half | LOOCV / unclipped | CLEAN | 74.898675 | 4.4162705 | [64.527176, 85.270175] | [62.541144, 87.256207] | CONF estimate lies inside the closed 95% prediction interval |
| P1c | sigma2_quarter | LOOCV / unclipped | CLEAN | 148.80711 | 6.2773662 | [134.06487, 163.54934] | [131.24189, 166.37232] | CONF estimate lies inside the closed 95% prediction interval |
| P1d | gamma | LOOCV / unclipped | CLEAN | 0.99043131 | 0.072483855 | [0.82020481, 1.1606578] | [0.78760825, 1.1932544] | CONF estimate lies inside the closed 95% prediction interval |
| P1e | CONF gamma 90% bootstrap CI contains 1 | LOOCV / unclipped | CLEAN | CI90 [0.87870395, 1.1146704] | N/A | N/A (categorical) | N/A (categorical) | CONF lo90 <= 1 <= hi90 |
| P2a | Delta_state = R_SB - ICC(C,1) | LOOCV / unclipped | PARTIALLY_KNOWN | -0.012066634 | 0.091228501 | [-0.22631446, 0.20218119] | [-0.26734064, 0.24320737] | CONF estimate lies inside the closed 95% prediction interval |
| P2b | S_counting | LOOCV / unclipped | PARTIALLY_KNOWN | 1.0209722 | 0.096997441 | [0.79317612, 1.2487682] | [0.7495556, 1.2923887] | CONF estimate lies inside the closed 95% prediction interval |
| P2c | CONF Delta_state 90% bootstrap CI contains 0 | LOOCV / unclipped | PARTIALLY_KNOWN | CI90 [-0.17493109, 0.12825751] | N/A | N/A (categorical) | N/A (categorical) | CONF lo90 <= 0 <= hi90 |
| P3 | Joint stored-position/adjacency absence criteria | LOOCV / clipped | CLEAN (ICC terms PARTIALLY_KNOWN) | All three criteria met on DEV | N/A | N/A (categorical) | N/A (categorical) | LMM p >= .01 AND ICC difference CI90 contains 0 AND lag1-minus-lag3 CI90 contains 0 |
| P3.p | Supporting P3 component: LMM_p | LOOCV / clipped | CLEAN | 0.22759967 | 0.21881442 | [-0.28628046, 0.74147981] | [-0.38468304, 0.83988239] | CONF estimate lies inside the closed 95% prediction interval |
| P3.ICCdiff | Supporting P3 component: ICC_difference | LOOCV / clipped | PARTIALLY_KNOWN | 0.00085248357 | 0.0030276776 | [-0.0062579406, 0.0079629077] | [-0.0076195112, 0.0093244783] | CONF estimate lies inside the closed 95% prediction interval |
| P3.lagdiff | Supporting P3 component: lag1_minus_lag3 | LOOCV / clipped | CLEAN | 0.050669865 | 0.10060262 | [-0.18559284, 0.28693257] | [-0.23083463, 0.33217436] | CONF estimate lies inside the closed 95% prediction interval |
| P4a | r_H_over_rstar | LOOCV / clipped | CLEAN | 0.63247087 | 0.04548715 | [0.52564545, 0.73929629] | [0.50518952, 0.75975222] | CONF estimate lies inside the closed 95% prediction interval |
| P4b | r_Q_over_rstar | LOOCV / clipped | CLEAN | 0.42713386 | 0.038560517 | [0.33657546, 0.51769225] | [0.31923449, 0.53503322] | CONF estimate lies inside the closed 95% prediction interval |
| P4c | r_H | LOOCV / clipped | CLEAN | 0.19710167 | 0.028181107 | [0.13091905, 0.26328429] | [0.11824579, 0.27595756] | CONF estimate lies inside the closed 95% prediction interval |
| P4d | r_Q | LOOCV / clipped | CLEAN | 0.13311095 | 0.020017746 | [0.086099776, 0.18012212] | [0.077097637, 0.18912426] | CONF estimate lies inside the closed 95% prediction interval |
| P5.1 | mean_contact_quality_vs_record_residual2 | LOOCV / clipped | CLEAN | -0.019640456 | 0.054596189 | [-0.14785822, 0.10857731] | [-0.17241056, 0.13312965] | CONF estimate lies inside the closed 95% prediction interval |
| P5.2 | mean_contact_quality_vs_half_difference2 | LOOCV / clipped | CLEAN | 0.016292276 | 0.036833349 | [-0.070209913, 0.10279447] | [-0.086774162, 0.11935871] | CONF estimate lies inside the closed 95% prediction interval |
| P5.3 | interpolated_fraction_vs_record_residual2 | LOOCV / clipped | CLEAN | -0.023497843 | 0.046425843 | [-0.13252776, 0.085532075] | [-0.15340583, 0.10641014] | CONF estimate lies inside the closed 95% prediction interval |
| P5.4 | interpolated_fraction_vs_half_difference2 | LOOCV / clipped | CLEAN | 0.028020832 | 0.043184023 | [-0.073395753, 0.12943742] | [-0.09281595, 0.14885761] | CONF estimate lies inside the closed 95% prediction interval |
| P5.5 | packet_counter_gaps_vs_record_residual2 | LOOCV / clipped | CLEAN | 0.0007733028 | 0.043888678 | [-0.10229815, 0.10384475] | [-0.12203523, 0.12358184] | CONF estimate lies inside the closed 95% prediction interval |
| P5.6 | packet_counter_gaps_vs_half_difference2 | LOOCV / clipped | CLEAN | -0.022243701 | 0.036893317 | [-0.10888672, 0.064399322] | [-0.12547794, 0.080990539] | CONF estimate lies inside the closed 95% prediction interval |
| P6.half.clipped | D_half | LOOCV / clipped | CLEAN | 0.017822402 | 0.0382462 | [-0.071997831, 0.10764264] | [-0.08919745, 0.12484225] | CONF estimate lies inside the closed 95% prediction interval |
| P6.quarter.clipped | D_quarter | LOOCV / clipped | CLEAN | 0.0089554403 | 0.016621612 | [-0.030079994, 0.047990875] | [-0.037554864, 0.055465745] | CONF estimate lies inside the closed 95% prediction interval |
| P6.label.clipped | A2 within-record tension label | LOOCV / clipped | CLEAN | NOT_DETECTED_WITHIN_RECORD | N/A | N/A (categorical) | N/A (categorical) | CONF A2 label equals registered DEV label |
| P7.D1.clipped | D1 | LOOCV / clipped | PARTIALLY_KNOWN | 0.069704105 | 0.023998376 | [0.013344526, 0.12606368] | [0.0025522658, 0.13685594] | CONF estimate lies inside the closed 95% prediction interval |
| P7.D4.clipped | D4 | LOOCV / clipped | PARTIALLY_KNOWN | 0.25762863 | 0.11196535 | [-0.0053191606, 0.52057643] | [-0.055670866, 0.57092813] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q09.clipped | rstar(q=0.9) / rstar(q=1); sensitivity | LOOCV / clipped | PARTIALLY_KNOWN | 0.92328619 | 0.092299661 | [0.70652277, 1.1400496] | [0.66501488, 1.1815575] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q09.clipped | r4(q=0.9) / r4(q=1); sensitivity | LOOCV / clipped | PARTIALLY_KNOWN | 0.96555658 | 0.087276494 | [0.76058994, 1.1705232] | [0.72134101, 1.2097721] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q08.clipped | rstar(q=0.8) / rstar(q=1); sensitivity | LOOCV / clipped | PARTIALLY_KNOWN | 0.76729265 | 0.10242659 | [0.52674639, 1.0078389] | [0.48068434, 1.053901] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q08.clipped | r4(q=0.8) / r4(q=1); sensitivity | LOOCV / clipped | PARTIALLY_KNOWN | 0.79305343 | 0.091974327 | [0.57705406, 1.0090528] | [0.53569247, 1.0504144] | CONF estimate lies inside the closed 95% prediction interval |
| P8.label.clipped | A4 counting-floor label; sensitivity | LOOCV / clipped | PARTIALLY_KNOWN | CONSISTENT_WITH_COUNTING | N/A | N/A (categorical) | N/A (categorical) | CONF A4 label equals registered DEV label |
| P6.half.unclipped | D_half | LOOCV / unclipped | CLEAN | 0.02908197 | 0.040255273 | [-0.06545652, 0.12362046] | [-0.083559635, 0.14172358] | CONF estimate lies inside the closed 95% prediction interval |
| P6.quarter.unclipped | D_quarter | LOOCV / unclipped | CLEAN | 0.01326397 | 0.017694944 | [-0.028292158, 0.054820097] | [-0.036249714, 0.062777654] | CONF estimate lies inside the closed 95% prediction interval |
| P6.label.unclipped | A2 within-record tension label | LOOCV / unclipped | CLEAN | NOT_DETECTED_WITHIN_RECORD | N/A | N/A (categorical) | N/A (categorical) | CONF A2 label equals registered DEV label |
| P7.D1.unclipped | D1 | LOOCV / unclipped | PARTIALLY_KNOWN | 0.068339945 | 0.023752414 | [0.012558, 0.12412189] | [0.0018763514, 0.13480354] | CONF estimate lies inside the closed 95% prediction interval |
| P7.D4.unclipped | D4 | LOOCV / unclipped | PARTIALLY_KNOWN | 0.24600445 | 0.10668269 | [-0.0045371561, 0.49654605] | [-0.052513208, 0.5445221] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q09.unclipped | rstar(q=0.9) / rstar(q=1); primary | LOOCV / unclipped | PARTIALLY_KNOWN | 0.91470623 | 0.087273777 | [0.70974597, 1.1196665] | [0.67049827, 1.1589142] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q09.unclipped | r4(q=0.9) / r4(q=1); primary | LOOCV / unclipped | PARTIALLY_KNOWN | 0.94203695 | 0.082539038 | [0.7481961, 1.1358778] | [0.71107764, 1.1729962] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q08.unclipped | rstar(q=0.8) / rstar(q=1); primary | LOOCV / unclipped | PARTIALLY_KNOWN | 0.79059746 | 0.097767433 | [0.56099312, 1.0202018] | [0.51702633, 1.0641686] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q08.unclipped | r4(q=0.8) / r4(q=1); primary | LOOCV / unclipped | PARTIALLY_KNOWN | 0.82116982 | 0.086538874 | [0.61793546, 1.0244042] | [0.57901825, 1.0633214] | CONF estimate lies inside the closed 95% prediction interval |
| P8.label.unclipped | A4 counting-floor label; primary | LOOCV / unclipped | PARTIALLY_KNOWN | CONSISTENT_WITH_COUNTING | N/A | N/A (categorical) | N/A (categorical) | CONF A4 label equals registered DEV label |
| P1a | r_HH | 46fold / unclipped | CLEAN | 0.0053090612 | 0.041717184 | [-0.09266269, 0.10328081] | [-0.11142324, 0.12204136] | CONF estimate lies inside the closed 95% prediction interval |
| P1b | sigma2_half | 46fold / unclipped | CLEAN | 77.315872 | 4.9798611 | [65.620794, 89.010949] | [63.381311, 91.250432] | CONF estimate lies inside the closed 95% prediction interval |
| P1c | sigma2_quarter | 46fold / unclipped | CLEAN | 149.18536 | 6.0838634 | [134.89756, 163.47316] | [132.1616, 166.20912] | CONF estimate lies inside the closed 95% prediction interval |
| P1d | gamma | 46fold / unclipped | CLEAN | 0.94826946 | 0.074417934 | [0.77350083, 1.1230381] | [0.74003449, 1.1565044] | CONF estimate lies inside the closed 95% prediction interval |
| P1e | CONF gamma 90% bootstrap CI contains 1 | 46fold / unclipped | CLEAN | CI90 [0.83434592, 1.0748189] | N/A | N/A (categorical) | N/A (categorical) | CONF lo90 <= 1 <= hi90 |
| P2a | Delta_state = R_SB - ICC(C,1) | 46fold / unclipped | PARTIALLY_KNOWN | -0.031997191 | 0.091078292 | [-0.24589226, 0.18189787] | [-0.28685088, 0.2228565] | CONF estimate lies inside the closed 95% prediction interval |
| P2b | S_counting | 46fold / unclipped | PARTIALLY_KNOWN | 1.0706076 | 0.098641104 | [0.83895146, 1.3022637] | [0.79459177, 1.3466234] | CONF estimate lies inside the closed 95% prediction interval |
| P2c | CONF Delta_state 90% bootstrap CI contains 0 | 46fold / unclipped | PARTIALLY_KNOWN | CI90 [-0.1890469, 0.11014425] | N/A | N/A (categorical) | N/A (categorical) | CONF lo90 <= 0 <= hi90 |
| P3 | Joint stored-position/adjacency absence criteria | 46fold / clipped | CLEAN (ICC terms PARTIALLY_KNOWN) | All three criteria met on DEV | N/A | N/A (categorical) | N/A (categorical) | LMM p >= .01 AND ICC difference CI90 contains 0 AND lag1-minus-lag3 CI90 contains 0 |
| P3.p | Supporting P3 component: LMM_p | 46fold / clipped | CLEAN | 0.26350723 | 0.23254181 | [-0.28261132, 0.80962579] | [-0.38718721, 0.91420168] | CONF estimate lies inside the closed 95% prediction interval |
| P3.ICCdiff | Supporting P3 component: ICC_difference | 46fold / clipped | PARTIALLY_KNOWN | 0.00063971289 | 0.0030145602 | [-0.0064399055, 0.0077193313] | [-0.0077955771, 0.0090750029] | CONF estimate lies inside the closed 95% prediction interval |
| P3.lagdiff | Supporting P3 component: lag1_minus_lag3 | 46fold / clipped | CLEAN | 0.065529986 | 0.098452029 | [-0.16568211, 0.29674209] | [-0.20995677, 0.34101674] | CONF estimate lies inside the closed 95% prediction interval |
| P4a | r_H_over_rstar | 46fold / clipped | CLEAN | 0.59989299 | 0.048702053 | [0.48551745, 0.71426853] | [0.46361575, 0.73617023] | CONF estimate lies inside the closed 95% prediction interval |
| P4b | r_Q_over_rstar | 46fold / clipped | CLEAN | 0.40605594 | 0.039153801 | [0.31410423, 0.49800765] | [0.29649645, 0.51561542] | CONF estimate lies inside the closed 95% prediction interval |
| P4c | r_H | 46fold / clipped | CLEAN | 0.18430252 | 0.028933515 | [0.1163529, 0.25225215] | [0.10334127, 0.26526378] | CONF estimate lies inside the closed 95% prediction interval |
| P4d | r_Q | 46fold / clipped | CLEAN | 0.12475081 | 0.019663922 | [0.078570583, 0.17093103] | [0.069727562, 0.17977405] | CONF estimate lies inside the closed 95% prediction interval |
| P5.1 | mean_contact_quality_vs_record_residual2 | 46fold / clipped | CLEAN | -0.028889521 | 0.05450004 | [-0.15688149, 0.099102445] | [-0.18139059, 0.12361154] | CONF estimate lies inside the closed 95% prediction interval |
| P5.2 | mean_contact_quality_vs_half_difference2 | 46fold / clipped | CLEAN | 0.044890774 | 0.035385521 | [-0.038211228, 0.12799278] | [-0.054124377, 0.14390592] | CONF estimate lies inside the closed 95% prediction interval |
| P5.3 | interpolated_fraction_vs_record_residual2 | 46fold / clipped | CLEAN | -0.016606734 | 0.0492646 | [-0.1323034, 0.099089933] | [-0.15445808, 0.12124461] | CONF estimate lies inside the closed 95% prediction interval |
| P5.4 | interpolated_fraction_vs_half_difference2 | 46fold / clipped | CLEAN | 0.0058060763 | 0.042734707 | [-0.094555301, 0.10616745] | [-0.11377344, 0.12538559] | CONF estimate lies inside the closed 95% prediction interval |
| P5.5 | packet_counter_gaps_vs_record_residual2 | 46fold / clipped | CLEAN | 0.015638362 | 0.044144231 | [-0.088033246, 0.11930997] | [-0.10788526, 0.13916198] | CONF estimate lies inside the closed 95% prediction interval |
| P5.6 | packet_counter_gaps_vs_half_difference2 | 46fold / clipped | CLEAN | -0.039571133 | 0.036348348 | [-0.12493431, 0.045792043] | [-0.14128045, 0.062138183] | CONF estimate lies inside the closed 95% prediction interval |
| P6.half.clipped | D_half | 46fold / clipped | CLEAN | 0.028506481 | 0.036678676 | [-0.057632462, 0.11464542] | [-0.074127153, 0.13114012] | CONF estimate lies inside the closed 95% prediction interval |
| P6.quarter.clipped | D_quarter | 46fold / clipped | CLEAN | 0.012553589 | 0.016310489 | [-0.025751183, 0.05085836] | [-0.033086139, 0.058193316] | CONF estimate lies inside the closed 95% prediction interval |
| P6.label.clipped | A2 within-record tension label | 46fold / clipped | CLEAN | NOT_DETECTED_WITHIN_RECORD | N/A | N/A (categorical) | N/A (categorical) | CONF A2 label equals registered DEV label |
| P7.D1.clipped | D1 | 46fold / clipped | PARTIALLY_KNOWN | 0.053481011 | 0.025221408 | [-0.0057508269, 0.11271285] | [-0.017093094, 0.12405512] | CONF estimate lies inside the closed 95% prediction interval |
| P7.D4.clipped | D4 | 46fold / clipped | PARTIALLY_KNOWN | 0.1905409 | 0.10821357 | [-0.063595922, 0.44467772] | [-0.11226042, 0.49334222] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q09.clipped | rstar(q=0.9) / rstar(q=1); sensitivity | 46fold / clipped | PARTIALLY_KNOWN | 0.91196978 | 0.09528907 | [0.68818581, 1.1357537] | [0.64533356, 1.178606] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q09.clipped | r4(q=0.9) / r4(q=1); sensitivity | 46fold / clipped | PARTIALLY_KNOWN | 0.96822263 | 0.09356318 | [0.74849187, 1.1879534] | [0.70641577, 1.2300295] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q08.clipped | rstar(q=0.8) / rstar(q=1); sensitivity | 46fold / clipped | PARTIALLY_KNOWN | 0.80045491 | 0.10964516 | [0.54295602, 1.0579538] | [0.49364773, 1.1072621] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q08.clipped | r4(q=0.8) / r4(q=1); sensitivity | 46fold / clipped | PARTIALLY_KNOWN | 0.84138561 | 0.097103142 | [0.61334134, 1.0694299] | [0.56967329, 1.1130979] | CONF estimate lies inside the closed 95% prediction interval |
| P8.label.clipped | A4 counting-floor label; sensitivity | 46fold / clipped | PARTIALLY_KNOWN | CONSISTENT_WITH_COUNTING | N/A | N/A (categorical) | N/A (categorical) | CONF A4 label equals registered DEV label |
| P6.half.unclipped | D_half | 46fold / unclipped | CLEAN | 0.030980619 | 0.039845995 | [-0.062596691, 0.12455793] | [-0.08051575, 0.14247699] | CONF estimate lies inside the closed 95% prediction interval |
| P6.quarter.unclipped | D_quarter | 46fold / unclipped | CLEAN | 0.0074782472 | 0.017115645 | [-0.032717412, 0.047673907] | [-0.040414454, 0.055370948] | CONF estimate lies inside the closed 95% prediction interval |
| P6.label.unclipped | A2 within-record tension label | 46fold / unclipped | CLEAN | NOT_DETECTED_WITHIN_RECORD | N/A | N/A (categorical) | N/A (categorical) | CONF A2 label equals registered DEV label |
| P7.D1.unclipped | D1 | 46fold / unclipped | PARTIALLY_KNOWN | 0.056332925 | 0.024759143 | [-0.001813297, 0.11447915] | [-0.01294768, 0.12561353] | CONF estimate lies inside the closed 95% prediction interval |
| P7.D4.unclipped | D4 | 46fold / unclipped | PARTIALLY_KNOWN | 0.19981924 | 0.10717448 | [-0.051877311, 0.45151578] | [-0.10007452, 0.49971299] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q09.unclipped | rstar(q=0.9) / rstar(q=1); primary | 46fold / unclipped | PARTIALLY_KNOWN | 0.91121039 | 0.0903895 | [0.69893294, 1.1234878] | [0.65828406, 1.1641367] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q09.unclipped | r4(q=0.9) / r4(q=1); primary | 46fold / unclipped | PARTIALLY_KNOWN | 0.95001419 | 0.08742695 | [0.74469422, 1.1553342] | [0.70537762, 1.1946508] | CONF estimate lies inside the closed 95% prediction interval |
| P8.rstar.q08.unclipped | rstar(q=0.8) / rstar(q=1); primary | 46fold / unclipped | PARTIALLY_KNOWN | 0.82770885 | 0.1055408 | [0.57984894, 1.0755688] | [0.5323864, 1.1230313] | CONF estimate lies inside the closed 95% prediction interval |
| P8.r4.q08.unclipped | r4(q=0.8) / r4(q=1); primary | 46fold / unclipped | PARTIALLY_KNOWN | 0.8653863 | 0.090944722 | [0.65180492, 1.0789677] | [0.61090636, 1.1198662] | CONF estimate lies inside the closed 95% prediction interval |
| P8.label.unclipped | A4 counting-floor label; primary | 46fold / unclipped | PARTIALLY_KNOWN | CONSISTENT_WITH_COUNTING | N/A | N/A (categorical) | N/A (categorical) | CONF A4 label equals registered DEV label |

## Label and verdict rules

**A2**: WITHIN_RECORD_NONINDEPENDENCE if the95% lower bound of D_half or D_quarter is>0; otherwise NOT_DETECTED_WITHIN_RECORD

**A4**: For each index separately: STEEPER_THAN_COUNTING if any95% upper bound of r(q)/r(1) is below q, across q=.9,.8 and rstar,r4; otherwise CONSISTENT_WITH_COUNTING. Unclipped is primary.

**A2_prespecified_interpretation**: A detected label is the prespecified diagnostic for within-record dependence in disjoint event subsets, placing tension within event detection/feature construction rather than requiring state differences between records. This rule alone does not establish a unique causal mechanism.

**counting_floor_assumption**: The supplied counting-floor model assumes the linear readout, expected features scaled by q, and counting noise unchanged in absolute units. Failure to reject this inequality does not establish the model.

**primary_readout**: LOOCV

**BUDGET_CONFIRMED**: At least3 of P1a-P1d pass AND P1e AND P2c AND P4a AND P4b pass

**BUDGET_NOT_CONFIRMED**: At most2 of P1a-P1d pass

**BUDGET_PARTIAL**: All other cases

**outside_verdict**: P3 and P5-P8; P2a and P2b also do not enter the specified verdict; 46-fold scored separately without changing primary verdict

## Frozen scorer, receipt and hashes

Run `src/stage2b_score.py --cohort DEV` for the dry run. Run with `--cohort CONF` only after registration. Set OPENBLAS_NUM_THREADS, OMP_NUM_THREADS, MKL_NUM_THREADS and VECLIB_MAXIMUM_THREADS to1 and PYTHONDONTWRITEBYTECODE=1; use the original venv and caffeinate -i. CONF bootstrap:2,000 resamples, seed20261025.

REGISTRATION_RECEIPT.txt must contain an https:// or http:// URL and ISO8601 timestamp with timezone; script refuses missing/invalid/future receipt. No receipt is created by Stage2a. Receipt SHA256 is currently **not available (file absent)**; Stage2b will record it. This document and JSON were not uploaded.

- `split_dev_conf.csv`: `69ece5f8828f7296374cc82125e0294ed8eef8cd5ff772f477cd2b0b8cb9606d`
- `configs/protocol_stage01.json`: `93b128b65b4a67431c0a29db603b33e58b87bd68e26bacb91f3b13371956ef3d`
- `configs/protocol_amendment_A1.json`: `ba4cc921cffe3d7a402fbda844929efe5daef828bda7a7f8d0a9fc89a2d5569e`
- `configs/protocol_amendment_A2.json`: `2bccbb52a77f8a4292d218f57694bfd8701a3018a0c240952d45e3dc266dafd0`
- `results/stage2_inputs.json`: `f374064865efcd53821d77a3cd76a4787b5bdfc9e2bae6225d17499bd712b2e7`
- `src/stage2b_score.py`: `a0958f212cb6e6d7f66912526ee8b5ebe4ee60214e583fd8ec18b999fba09154`
- `src/common.py`: `9d614cb5ec4ef4d04ce763162c0c22533521d0c86974edb900008d262613043e`
- `src/inference.py`: `235f4f3af672efa40ec914e45baaf8ee55863ab9651fde6c290ce20d26ec002f`
- `src/analysis.py`: `a7c6987d8147089f0f11d527e5b1ef861832ddc239bc7774474e0a8205f29c35`
- `src/a2_features.py`: `088180188bc5c6f6af78940b0534e1e49ae340db52e85e8fc54aea2c97057546`
- `src/stage1_math.py`: `f32f51a4ee8d242aac751544dd3156b09abf4439aeb9611bb730467d1bdecadc`
- `V6 protocol_stage2a.json`: `1e247513581b6ff25117b48c010128d0815dc78656e3f79fc32c5f7c41ec0248`
- `DEV bootstrap`: `09de7a6285e10964f32c86303f200cba460911e6698cd2e32008c2833879dd81`

## Limitations

- DEV is half of a heavily explored cohort; one visit and one SDS per person.
- Recording order is unverified; no confirmed temporal or state interpretation.
- Intervals resample people only; one fixed thinning realization per record and q is conditioned upon.
- No multiple-testing adjustment or equivalence test is added; failure to detect is not proof of zero.
- Full rstar is position-centred; pooled subrecord accuracy is uncentred, as fixed before computations.
- The scorer checks registration bytes and hashes of every executable dependency; software/environment versions are recorded separately.

Timestamp (UTC): 2026-10-24T11:50:39Z
