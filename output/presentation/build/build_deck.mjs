import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import sharp from 'sharp';
import { Presentation, PresentationFile } from '@oai/artifact-tool';

const root = 'C:/Users/naman/Desktop/planetx-spherex-lab';
const build = `${root}/output/presentation/build`;
const output = `${root}/output/presentation/deliverables`;
const skill = 'C:/Users/naman/.codex/plugins/cache/openai-primary-runtime/presentations/26.1112.11217/skills/presentations';
const runtime = 'C:/Program Files/WindowsApps/OpenAI.CodexPrimaryRuntime.v26-1112-737-0_26.1112.737.0_x64__3k8sg7r9htsxt/dependencies';
process.env.RUNTIME_NODE_MODULES = `${runtime}/node/node_modules`;
process.env.RUNTIME_NODE = `${runtime}/node/bin/node.exe`;
process.env.RUNTIME_PYTHON = `${runtime}/python/python.exe`;
process.env.PYTHONUTF8 = '1';
const font = 'Segoe UI';
const C = { bg:'#0A1118', white:'#F5F5F0', muted:'#B4C1CC', orange:'#F5A14C', blue:'#75BFF3', paper:'#F5F3ED', ink:'#10212C', grey:'#4E626E' };
await fs.mkdir(`${build}/assets`, {recursive:true});
await fs.mkdir(`${build}/previews`, {recursive:true});
await fs.mkdir(output, {recursive:true});

async function crop(file, name, box) {
  const out = `${build}/assets/${name}.png`;
  await sharp(`${root}/${file}`).extract(box).png().toFile(out);
  return out;
}
const assets = {
  cover: await crop('output/playwright/current-eps-pair.png','cover-field',{left:158,top:640,width:420,height:420}),
  spectral: await crop('output/playwright/current-spectral-full.png','spectral-match',{left:269,top:574,width:635,height:638}),
  viewer: await crop('output/playwright/current-comparison.png','comparison-ui',{left:120,top:88,width:1200,height:912}),
  significance: await crop('artifacts/epoch_lab/demos/barnard-pair/comparison.png','barnard-significance',{left:1408,top:64,width:445,height:400}),
  series: `${root}/artifacts/epoch_lab/demos/eps-indi-b-series/series.png`,
};

const p = Presentation.create({slideSize:{width:1600,height:900}});
function text(s, name, value, x,y,w,h,size=30,color=C.white,bold=false, extra={}) {
  const shape = s.shapes.add({name,geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  shape.text = value;
  shape.text.style = {typeface:font,fontSize:size,color,bold,autoFit:'none',wrap:'square',insets:0,verticalAlignment:'top',...extra};
  return shape;
}
async function img(s, source, alt, x,y,w,h,fit='contain') {
  const bytes = await fs.readFile(source);
  s.images.add({blob:bytes,contentType:'image/png',alt,fit,position:{left:x,top:y,width:w,height:h}});
}
function base(n, light=false) {
  const s=p.slides.add(); s.background.fill=light?C.paper:C.bg;
  text(s,'slide-number',String(n).padStart(2,'0'),1490,839,40,30,19,light?C.grey:C.muted);
  return s;
}
function title(s,t,light=false){text(s,'slide-title',t,80,65,1430,90,58,light?C.ink:C.white,true);}
function notes(s,body,sources){s.speakerNotes.text=`${body}\n\nSources\n${sources.join('\n')}`;}

// 1 — Cover, using a real SPHEREx field rather than an illustrative sky.
{
 const s=base(1);
 await img(s,assets.cover,'Real SPHEREx DR1 image of the epsilon Indi B field, cropped from the recorded one-year pair viewer',872,0,728,900,'cover');
 text(s,'event','NASA SPACE APPS CHALLENGE 2026',80,74,750,40,25,C.muted);
 text(s,'project-name','SPHEREx\nEpoch Lab',76,242,790,255,104,C.white,true);
 text(s,'project-purpose','Compare the sky across time\nat the same wavelength',80,550,750,110,39,C.white);
 text(s,'challenge','Planet X and SPHEREx',80,725,750,48,30,C.orange);
 text(s,'team-name','Kugelblitz Vanguard',80,790,750,36,25,C.white,true);
 text(s,'team-member','Team member: Naman P Menasagi',80,834,750,32,22,C.muted);
 notes(s,
 'Team Kugelblitz Vanguard: Naman P Menasagi. SPHEREx Epoch Lab is a working browser interface for inspecting real SPHEREx observations across time. Its central scientific choice is to compare only measurements that sampled the same wavelength. The project includes archive screening, an interactive image viewer, validation reports and multi-epoch light curves and motion analysis. The image is an actual processed DR1 field from the recorded epsilon Indi B one-year pair. The bright flagged core visible in this field is excluded by the pipeline mask. The challenge name does not imply a new-object discovery. No discovery is claimed.',
 ['Repository: README.md; docs/SUBMISSION.md; static/index.html; static/app.js',
  'Image: current app, #runs/eps-indi-b-year-pair; underlying real archive files in artifacts/epoch_lab/demos/eps-indi-b-year-pair/',
  'Team name and member supplied by user for this deck',
  'Official challenge description, About the Challenge section: https://www.spaceappschallenge.org/2026/find-a-team/lsp26/ (accessed 2026-11-15)',
  'NASA Space Apps 2026 event: https://www.spaceappschallenge.org/']);
}

// 2 — Scientific problem and the actual restricted-overlap example.
{
 const s=base(2); title(s,'The comparison trap');
 text(s,'challenge-brief','The challenge asks for a public web tool\nto inspect changes in SPHEREx sky images.',80,191,790,100,34,C.white);
 text(s,'science-mechanism','SPHEREx samples 102 infrared colours.\nA star’s detector position changes the\nwavelength measured on each visit.',80,326,790,145,34,C.muted);
 text(s,'false-variability','Subtracting different colours can create\nfalse variability.',80,495,780,95,34,C.white,true);
 text(s,'matched-percent','7.5%',80,630,240,106,88,C.orange,true);
 text(s,'matched-explanation','of overlapping pixels match in this\nRW Aurigae pair: 1,483 of 19,710.',330,653,535,85,30,C.white);
 text(s,'spectral-policy','Lab policy: |Δλ| ≤ 0.1 × bandwidth',80,793,780,45,27,C.muted);
 await img(s,assets.spectral,'RW Aurigae spectral-match map: the narrow blue strip passes the wavelength gate while orange areas differ in wavelength',920,190,590,592);
 text(s,'map-caption','Blue passes the wavelength gate.\nOrange samples different colours.',920,793,560,62,24,C.muted);
 notes(s,
 'The challenge is to make temporal changes in SPHEREx images easy for the public to inspect. SPHEREx uses linear variable filters: the band center varies across each detector. Seeing the same sky twice therefore does not guarantee comparable wavelengths. The lab reads CWAVE and CBAND calibration products and gates every compared pixel. In the recorded RW Aurigae comparison, observations on 2025-09-03 and 2026-03-04 overlap on 19,710 valid pixels, but only 1,483 satisfy the spectral policy, or 7.524%. The displayed map is the current application’s spectral-match view of that real pair. Its blue strip identifies the spectrally comparable area. The 10% bandwidth tolerance is the lab’s policy, not a NASA calibration specification. This pair is PARTIAL because its noise gates could not both be verified, and its target core is flagged.',
 ['NASA/JPL filter design and 102 wavelengths: https://www.jpl.nasa.gov/press-kits/spherex/mission-overview/spacecraft/',
  'Official CWAVE/CBAND guidance: https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/',
  'Official challenge wording: https://www.spaceappschallenge.org/2026/find-a-team/lsp26/',
  'Project evidence: artifacts/epoch_lab/demos/rw-aurigae-strip/report.json, spectral; difference_imaging.py, POLICY and compare; docs/METHODOLOGY.md',
  'Image: current app #runs/rw-aurigae-strip, Spectral match view, captured 2026-11-15']);
}

// 3 — Public-facing workflow with a fresh product screenshot.
{
 const s=base(3); title(s,'A web workflow for temporal exploration');
 const steps=[
  ['01','Choose a field','Preset, object name or RA/Dec.'],
  ['02','Find matching epochs','Rank pairs by matched sky area.'],
  ['03','Compare images','Blink, overlay and significance.'],
  ['04','Inspect the evidence','Residuals, validation and downloads.'],
 ];
 for(let i=0;i<steps.length;i++){
  const y=222+i*143;
  text(s,`step-${i}-number`,steps[i][0],80,y,75,38,25,C.orange);
  text(s,`step-${i}-heading`,steps[i][1],80,y+36,425,49,32,C.white,true);
  text(s,`step-${i}-detail`,steps[i][2],80,y+88,425,70,26,C.muted);
 }
 await img(s,assets.viewer,'Current Epoch Lab recorded Barnard field comparison with validation badges and side-by-side image viewer',585,180,935,640);
 notes(s,
 'A user selects a preset, resolves an object name or enters coordinates. Screening maps archive footprints to detector coordinates and nominates same-detector pairs and wavelength groups. Comparisons download real cutouts and re-check each observation’s own calibration. The image viewer offers side-by-side, overlay, blink, difference, significance and spectral-match views, with zoom, pan and pixel readout. Results expose residual classifications, validation gates, per-star fluxes and provenance files. Groups of three or more candidate epochs can launch a time series. The screenshot is a recorded real-data run, clearly labelled by the app. Recorded runs open offline. New searches and analyses require the local Python server and available upstream services. The live route is implemented; this screenshot does not assert that IRSA is currently reachable.',
 ['Implementation: static/index.html; static/app.js; epoch_lab_api.py; spectral_matching.py; epoch_lab.py',
  'Routes: POST /api/lab/screen, /api/lab/compare, /api/lab/series; GET /api/lab/demos',
  'Image: current local app #runs/barnard-pair, captured 2026-11-15',
  'Workflow and availability: docs/ARCHITECTURE.md; build_static_demo.py']);
}

// 4 — Native editable evidence table, including real integration boundaries.
{
 const s=base(4); title(s,'Real data, traceable results');
 text(s,'data-layer-summary','Python science pipeline and FastAPI, with a plain JavaScript interface',80,163,1400,50,31,C.muted);
 const values=[
  ['Source','Implemented use','Availability and limits'],
  ['NASA/IPAC IRSA','SPHEREx DR1 search and FITS cutouts','Live requests or validated local caches'],
  ['IRSA DataLink','Per-observation CWAVE and CBAND maps','Missing calibration fails that epoch'],
  ['ESA Gaia DR3','Reference stars and residual associations','Stored queries or live TAP. Absence is labelled'],
  ['JPL Horizons','SPHEREx observer vectors for parallax factors','Geometric vectors with response provenance'],
 ];
 const table=s.tables.add({rows:5,columns:3,left:80,top:255,width:1440,height:410,columnWidths:[305,565,570],values});
 table.styleOptions={headerRow:false,bandedRows:false};
 table.borders.assign({style:'solid',fill:'#34424C',width:1});
 for(let r=0;r<5;r++){
  table.rows[r].height=r===0?62:87;
  for(let c=0;c<3;c++){
   const cell=table.getCell(r,c); cell.fill=C.bg;
   cell.text.style={typeface:font,fontSize:r===0?25:27,color:r===0?C.orange:C.white,bold:r===0||c===0,insets:{left:14,right:15,top:13,bottom:10},verticalAlignment:'middle'};
  }
 }
 text(s,'provenance','Each result records its science files, calibration files and code hashes.',80,708,1400,55,32,C.white,true);
 text(s,'outputs','FITS, JSON, CSV and PNG downloads keep the evidence inspectable.',80,776,1380,60,30,C.muted);
 notes(s,
 'The application uses real NASA/IPAC IRSA SPHEREx DR1 Level 2 images and observation-associated spectral WCS calibration. Screening uses a footprint model to estimate the detector location, while image comparison re-reads each cutout’s own DataLink association. Processing masks flagged pixels, fits backgrounds, cross-convolves measured effective PSFs, registers images and applies the per-pixel spectral gate. Gaia DR3 reference queries are reused only if their sky circle covers the field. If Gaia is unavailable the report records that absence and limits source classification; it does not fabricate a catalog. Horizons vectors supply the observer positions for annual-parallax factors. Direct science/calibration UUID equality remains recorded as unverified: DataLink association and file hashes are the available provenance. Missing calibrations are never replaced with another epoch’s files. Time series can continue as PARTIAL when enough other epochs remain. Synthetic skies are used for tests and artificial injections are explicitly labelled.',
 ['Project implementation: spherex.py; calibration_resolver.py; fits_pipeline.py; spectral_matching.py; difference_imaging.py; epoch_lab.py; epoch_series.py; reference_catalog.py; spacecraft_ephemeris.py',
  'Architecture and integration inventory: docs/ARCHITECTURE.md; docs/DATA_SOURCES.md',
  'IRSA access: https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/',
  'IRSA spectral products: https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/',
  'Horizons endpoint used in code: https://ssd.jpl.nasa.gov/api/horizons.api']);
}

// 5 — Real comparison with injection evidence clearly distinguished.
{
 const s=base(5); title(s,'A validated comparison of the Barnard field');
 text(s,'baseline','31 days apart at 0.83–0.86 µm',80,165,1000,48,32,C.muted);
 text(s,'gate-number','5 of 5',80,249,750,100,82,C.orange,true);
 text(s,'gate-label','validation gates passed',80,348,750,58,38,C.white,true);
 text(s,'sample-size','23,265 matched pixels\n909 stars calibrate the source noise model',80,433,750,105,32,C.muted);
 text(s,'injection-result','27 of 31 injected sources recovered',80,577,760,60,36,C.white,true);
 text(s,'injection-disclosure','Synthetic S/N-10 sources added to the real\ndifference image. Detection threshold: 5σ.',80,640,760,85,28,C.muted);
 text(s,'science-result','No vetted change detected in this field',80,754,1340,58,37,C.white,true);
 text(s,'saturated-target','Barnard’s Star itself is saturated and excluded.',80,825,1350,40,25,C.muted);
 await img(s,assets.significance,'Real Barnard field significance image from the recorded validated comparison; units are sigma and orange contours mark the compared region',925,240,590,530);
 notes(s,
 'The saved Barnard field pair compares 2025W35_2A_0474_1/D1 with 2025W40_1A_0441_1/D1 over 31.0497 days. All five implemented status gates pass: spectral overlap, photometric consistency, blank-sky pixel noise, source noise model and registration. The significance figure is a cropped original pipeline output, preserving its scale and contour information. Median registration offsets are 0.036 and 0.030 arcsec using 631 stars. The fitted photometric scale B/A is 1.0033 ± 0.0037, and no rescaling was needed. Blank-sky spread before uncertainty inflation was 1.369; source z-score robust spread is 0.999 after empirical calibration. Five of 909 measured stars exceeded 5 sigma, and all five failed vetting due to blends or edges. There are zero vetted candidates. Artificial source injection recovered 27 of 31 available S/N-10 trials. That is a diagnostic on this real field, not a survey-wide completeness measurement and not validation on a real variable star. Barnard’s Star’s own core is masked by quality flags, so the result concerns the surrounding field.',
 ['Saved evidence: artifacts/epoch_lab/demos/barnard-pair/report.json; comparison.fits; comparison.png; candidates.csv',
  'Independent arithmetic: verify_results.py; docs/VERIFICATION.md, Barnard field section',
  'Processing and gate logic: difference_imaging.py, compare, _injection_recovery, _vet_variable',
  'Image: artifacts/epoch_lab/demos/barnard-pair/comparison.png, significance panel']);
}

// 6 — Original scientific output, with its actual scope and caveats.
{
 const s=base(6,true); title(s,'ε Indi B: a known mover recovered',true);
 text(s,'epoch-baseline','Five epochs at 1.107 µm over one year',80,163,1430,52,33,C.grey);
 await img(s,assets.series,'Original SPHEREx epsilon Indi B output: same-wavelength light curve using three accepted photometric epochs and trajectory measured from five astrometric epochs with error bars and catalog prediction',80,238,1440,528);
 text(s,'motion-significance','11.3σ motion detection',80,785,710,58,43,C.ink,true);
 text(s,'catalog-agreement','1.1σ from catalog proper motion',815,785,710,58,37,C.ink,true);
 text(s,'scientific-limits','Flux is consistent with constant in 3 usable epochs. Parallax remains poorly constrained.',80,852,1400,36,25,C.grey);
 notes(s,
 'The recorded epsilon Indi B series contains five astrometric epochs from 2025-04-24 through 2026-04-24 at about 1.107 microns. The proper-motion-plus-parallax fit measures east and north motion of (4.46 ± 0.45, −2.33 ± 0.43) arcsec per year. The covariance-based proper-motion significance is 11.309 sigma, and the difference from the independent SIMBAD preset proper motion is 1.116 sigma. This is verification on one known mover. The plotted catalog trajectory uses JPL Horizons SPHEREx observer vectors. The trajectory residuals exceed the nominal model noise: the fit applies a Birge factor of 3.402. Parallax is 0.611 ± 0.164 arcsec versus catalog 0.271 arcsec and is not reliably constrained. Only three epochs have enough wavelength-matched field stars for a photometric zero point. Those yield a weighted constant flux of 14.488 mJy with p=0.484. The other two photometric epochs are explicitly excluded. The figure is the original saved scientific product, with axes, uncertainties and legend preserved.',
 ['Source data and statistics: artifacts/epoch_lab/demos/eps-indi-b-series/series.json, astrometry.astrometric_solution and photometry.constant_flux_test',
  'Original image: artifacts/epoch_lab/demos/eps-indi-b-series/series.png',
  'Independent recomputation: verify_results.py; docs/VERIFICATION.md, epsilon Indi B series section',
  'Implementation: epoch_series.py; spacecraft_ephemeris.py; targets.py (SIMBAD links and preset values)']);
}

// 7 — Honest readiness, intended benefit and a concrete next milestone.
{
 const s=base(7); title(s,'Ready for exploration, with explicit limits');
 text(s,'today-heading','Works today',80,217,680,54,35,C.orange,true);
 text(s,'today-details','Recorded analyses open offline.\nLocal live analysis connects to\nthe NASA/IRSA archive.',80,296,690,160,36,C.white);
 text(s,'verification-count','81 passed\n0 failed\n3 unverified',80,505,675,195,48,C.white,true);
 text(s,'verification-scope','84 independent checks of the saved products',80,720,710,50,27,C.muted);
 text(s,'next-heading','Next scientific milestone',890,217,630,54,35,C.orange,true);
 text(s,'next-step','Verify sensitivity on an\nunsaturated variable star.',890,296,630,150,40,C.white,true);
 text(s,'limits','Archive outages, saturation and\nsparse wavelength overlap still\nlimit which fields can be analysed.',890,492,630,155,32,C.muted);
 text(s,'discovery-boundary','No discovery is claimed. Candidates need independent follow-up.',80,790,1400,47,31,C.white);
 const link=text(s,'repository','Source code and run instructions: github.com/Namanpm/planetx-spherex-lab',80,850,1390,35,25,C.blue);
 link.text.get('github.com/Namanpm/planetx-spherex-lab').link={uri:'https://github.com/Namanpm/planetx-spherex-lab',isExternal:true};
 notes(s,
 'The project makes technically demanding archive comparisons explorable through a public-facing interface, while keeping validity and provenance visible. The local server supports live screening, comparisons and time series when upstream archives are available. Recorded products use real retrieved archive data and require no network. The static builder provides a server-free copy of the recorded runs and the repository contains a GitHub Pages publishing workflow; no deployed public-site URL is asserted here. A fresh read-only invocation of verify_results.py --no-markdown recomputed 84 checks on the saved products: 81 PASS, 0 FAIL and 3 UNVERIFIED. These are checks of this recorded evidence, not 84 independent discoveries or a broad accuracy benchmark. Motion has been checked on one known mover. Real variable-star sensitivity is not yet demonstrated because the variable presets saturate and archive calibration requests failed during attempted runs. The next milestone is a completed, calibrated series of an unsaturated variable with enough matching epochs. The source-code link comes from the repository’s configured Git remote.',
 ['Readiness and limits: README.md; docs/SUBMISSION.md; docs/METHODOLOGY.md; docs/DATA_SOURCES.md',
  'Recorded-product verification: verify_results.py; output/presentation/build/verification.json, recomputed 2026-11-15',
  'Static export and publication workflow: build_static_demo.py; .github/workflows/pages.yml',
  'Repository: https://github.com/Namanpm/planetx-spherex-lab']);
}

await fs.writeFile(`${build}/presentation.json`,JSON.stringify(p.toProto()));
for(const [i,s] of p.slides.items.entries()) {
 const name=`slide-${String(i+1).padStart(2,'0')}`;
 console.log(`Rendering ${name}`);
 const blob=await p.export({slide:s,format:'png',scale:1});
 await fs.writeFile(`${build}/previews/${name}.png`,new Uint8Array(await blob.arrayBuffer()));
 const layout=await s.export({format:'layout'});
 await fs.writeFile(`${build}/previews/${name}.json`,await layout.text());
}
const candidate=`${build}/candidate.pptx`;
await (await PresentationFile.exportPptx(p)).save(candidate);
console.log('Draft PPTX exported');
const {finalizePresentation}=await import(pathToFileURL(`${skill}/container_tools/artifact_tool_utils.mjs`).href);
const result=await finalizePresentation({
 workspaceDir:root,candidatePath:candidate,finalPath:`${output}/SPHEREx_Epoch_Lab_Space_Apps_2026_UPDATED.pptx`,
 pythonExecutable:`${runtime}/python/python.exe`,
 integrityValidatorPath:`${skill}/container_tools/inspect_presentation_package_integrity.py`,
 layoutValidatorPath:`${skill}/container_tools/inspect_presentation_layout_geometry.py`,
 layoutArgs:['--expected-slide-size-emu','15240000,8572500','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','4'],
 explicitTotalSlideCount:7,requiredNativeTableOwnerSlides:[4],requiredNativeChartOwnerSlides:[],
 fontPolicy:{basis:'design',families:[font]},verifyArtifactToolImport:true,
 receiptPath:`${build}/validation-updated.json`,
});
console.log(JSON.stringify(result));
try {
 const pdf=await p.export({format:'pdf'});
 await fs.writeFile(`${output}/SPHEREx_Epoch_Lab_Space_Apps_2026_UPDATED.pdf`,new Uint8Array(await pdf.arrayBuffer()));
 console.log('PDF exported');
} catch(error) {console.log(`Optional PDF not exported: ${error.message}`);}
