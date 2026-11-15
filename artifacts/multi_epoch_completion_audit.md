# Multi-epoch completion audit

Specification: user attachment `5fd178dd-30da-4a54-9ab1-e71f3cb2b233/pasted-text-1.txt`.

| Requirement | Current implementation | Verification / boundary |
|---|---|---|
| Inspect/reuse existing system | Existing native fits, cache and empirical pair tests retained | Baseline137 tests before changes; final199 passing |
| Select2..8, 3modes | observation_selection.py, multi_epoch.py, API/UI | Tests2/5/8=>1/10/28; explicitIDs never substituted |
| Scientific rules | Existing pair policy; timeline independent applicability gate; explicit failures | Original520measurements,26pairrecords,confirmationexact; no fakeepochs |
| Readiness/live/replay | Direct selectedFITS inspection, calibrationchecksum/mapvalidation, inputhashes; boundedretry | Orchestration testswrongdetector/missingcal/noHTTPfallback; local502/503retrytests |
| Prepareonce/caches | run_selected uses preparedobservationcache, calibrationLRU, batchedareas; replay inspectedonce/selectedfile | Freshfit2/5/8 benchmark; existingcachetests retained |
| Sourceidentity | Catalog logical IDs, PMpredictions, availableuncertainties, explicitforcedpositionsemantics | Proper-motion regression; mismatchedforcedposition rejection; uncertainty spatialcandidate tests |
| Measureonce | Observation packets shared acrossallpairs; explicitmeasuredsnapshot replay | Freshall6 benchmark runs match savedmeasurement/rejectiondictsexactly |
| Allpairs/bestpair | Every mathematicalpair recorded; sameepoch/crossdetector rejection; existingvalidationrankingandoneconfirmation | 28records include26unchangedoriginalrecords; selectedpairconfirmationfailure stopswithoutfallback |
| Jointtimeseries | Simultaneoussource-nuisance/epochintercept/colorWLS; fullsharedcov; applicabilitygate | Syntheticnull/injection/rank/sharedcov/missingepochs tests; real5INSUFFICIENT |
| Independentreferencepartitions | Existingglobalhashgroups; jointglobalstampdisjointness; target excluded | Testsconfirmationdoesnotchangefit/subset; stampoverlap tests |
| Scientificoutputs | readiness/source/pairFITS, actualdateplots, conditionaltimelinecovFITS, manifest | EveryrealreplayFITSchecksumandartifacthash passes; syntheticvalidtimelineartifacttest |
| Backend | Boundedbackgroundworker, persistedjobstatus, input/pathguards, results/readiness/pairs/timeline/artifacts | APIinput/schema/busy/path/failuretests; actualbrowserjobcompleted |
| Frontend | ServedHTML/JS, target/coordinates, select2..8 cards, modes, actualstagepolling, outputs/rejections, reducedmotion | Chromeactualselection5=>10pairs; IABlive226records; HTMLdropdownregression |
| Performance | Local/networkphasecounters, peakRSS, actualfreshfitcoldprocess/warmarrays +snapshotruns | validation_v1 timings; discoveryexcludedexplicitly;286.49MiBcumulativeRSS |
| Tests | Existing137 +62new tests |199pass;2preexistingwarnings; noassertionsweakened |
| Realvalidation | Baselineexact regression; real5selection; fresh2/5/8 photometry | validation_final/validation.json true +validation_v1/validation.json true |
| Integratedimplementation | Reusablefunctions, CLI, API, UI, artifacts | ActualcommandsinREADME; modulescompile; diffwhitespacecheckpasses |
| Acceptancecriteria | Allengineeringmodesexposed, validrejectionsvisible, legacycompatible | No acceptedreal5timelineclaim; scientificinsufficiencyisexpectedvalidoutput |
| Delivery | READMEcommands, exactIDs/results/timings/limitations | Finalreportlinksvalidationrecords and sourcefiles |

Scientific boundaries: accepted pair counts mean reference validation consistency, not independent confirmation of every pair or physical response equality. Real common-scale applicability is unverified; no RW Aurigae variability or moving-object detection is claimed. Measured astrometric trajectories are unsupported because current photometry is catalog-forced; positional matching returns only candidates. Unknown systematics are not manufactured. Snapshot hashes describe current replay inputs, not proof of an unavailable historical science checksum. Cold-process benchmarks reuse raw disk products; they are not full cold-download measurements. Memory is a cumulative process high-water mark.
