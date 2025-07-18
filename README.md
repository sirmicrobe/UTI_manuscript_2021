# Performance of Conventional Urine Culture Compared to 16S rRNA Gene Amplicon Sequencing in Children with Suspected Urinary Tract Infection
The published manuscript is here:
https://journals.asm.org/doi/full/10.1128/spectrum.01861-21 <br>
Abstract:<br>
<p>Because some organisms causing urinary tract infection (UTI) may be difficult to culture, examination of bacterial gene sequences in the urine may provide a more accurate view of bacteria present during a UTI. Our objective was to estimate how often access to 16S rRNA gene amplicon sequencing alters diagnosis and/or clinical management. The study was designed as a cross-sectional study of a convenience sample of children with suspected UTI. The setting was the emergency department or outpatient clinic at six pediatric centers. Participants included children 2 months to 10 years of age suspected of UTI. We categorized the results of urine culture as follows: “likely UTI” (≥100,000 CFU/ml of a single uropathogen), “possible UTI” (10,000 to 99,000 CFU/ml of a uropathogen or ≥100,000 CFU/ml of a single uropathogen plus other growth), and “unlikely UTI” (no growth or growth of nonuropathogens). Similarly, we categorized the results of 16S rRNA gene sequencing into the same three categories using the following criteria: likely UTI (≥90% relative abundance of a uropathogen), possible UTI (50 to 89% relative abundance of a uropathogen), and unlikely UTI (remainder of samples). The main study outcome was concordance between conventional culture results and 16S rRNA gene sequencing. Concordance between the two methods was high in children with likely and unlikely UTI by conventional culture (95% and 87%, respectively). In children with possible UTI according to conventional culture, 71% had a single uropathogen at a relative abundance of ≥90% according to 16S rRNA gene sequencing data. Concordance between conventional culture and 16S rRNA gene amplicon sequencing appears to be high. In children with equivocal culture results, 16S rRNA gene results may provide information that may help clarify the diagnosis.</p> <br>
<br>




This repository contains R code for 16S rRNA gene amplicon sequence analysis on 118 kids with suspected UTI. <br>
<br>










<p>Amplicon sequencing methods from Earth Microbiome Project<br>
https://www.protocols.io/view/emp-16s-illumina-amplicon-protocol-kqdg3dzzl25z/v2<br>

From: https://earthmicrobiome.org/protocols-and-standards/16s/<br>
16S V4 amplification primers<br>
The current primers have been modified from the original 515F–806R primer pair (Caporaso et al., 2011) in the following ways:<br>
<br>
Barcodes are now on the forward primer 515F (Parada et al., 2016). This enables the usage of various reverse primer constructs to obtain longer amplicons, for example the V4–V5 region using reverse primer 926R (Quince et al., 2011; Parada et al., 2016).
Degeneracy was added to both the forward and reverse primers to remove known biases against Crenarachaeota/Thaumarchaeota (515F, also called 515F-Y, Parada et al., 2016) and the marine and freshwater Alphaproteobacterial clade SAR11 (806R, Apprill et al., 2015).
The primer sequences without linker, pad, barcode, or adapter are as follows:<br>

Updated sequences: 515F (Parada)–806R (Apprill), forward-barcoded:<br>
FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT<br>
Original sequences: 515F (Caporaso)–806R (Caporaso), reverse-barcoded:<br>
FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT <p> <br>
