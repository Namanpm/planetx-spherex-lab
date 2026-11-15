import fs from 'node:fs/promises';
import {PresentationFile, FileBlob} from '@oai/artifact-tool';
import sharp from 'sharp';
const base='C:/Users/naman/Desktop/planetx-spherex-lab/output/presentation';
const source=`${base}/deliverables/SPHEREx_Epoch_Lab_Space_Apps_2026.pptx`;
const out=`${base}/build/visual-qa`;
await fs.mkdir(out,{recursive:true});
const p=await PresentationFile.importPptx(await FileBlob.load(source));
for (const [i,s] of p.slides.items.entries()) {
 const blob=await p.export({slide:s,format:'png',scale:1});
 await fs.writeFile(`${out}/slide-${i+1}.png`,new Uint8Array(await blob.arrayBuffer()));
 console.log(`Rendered slide ${i+1}`);
}
const tiles=[];
for(let i=0;i<p.slides.items.length;i++){
 const input=await sharp(`${out}/slide-${i+1}.png`).resize(640,360).png().toBuffer();
 tiles.push({input,left:(i%2)*660+20,top:Math.floor(i/2)*380+20});
}
await sharp({create:{width:1340,height:1540,channels:3,background:'#25313b'}}).composite(tiles).png().toFile(`${out}/seven-slide-overview.png`);
console.log(`All ${p.slides.items.length} slides rendered from the delivered PPTX`);
