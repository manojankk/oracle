MERGE INTO bonuses D
   USING (SELECT employee_id, salary, department_id FROM employees
   WHERE department_id = 80) S
   ON (D.employee_id = S.employee_id)
   WHEN MATCHED THEN UPDATE SET D.bonus = D.bonus + S.salary*.01
     DELETE WHERE (S.salary > 8000)
   WHEN NOT MATCHED THEN INSERT (D.employee_id, D.bonus)
     VALUES (S.employee_id, S.salary*.01)
     WHERE (S.salary <= 8000);
	 
	 
 MERGE INTO dest_tab a
    USING source_tab b
    ON (a.object_id = b.object_id)
    WHEN MATCHED THEN
      UPDATE SET
        owner       = b.owner,
        object_name = b.object_name,
        object_type = b.object_type
    WHEN NOT MATCHED THEN
      INSERT (object_id, owner, object_name, object_type)
      VALUES (b.object_id, b.owner, b.object_name, b.object_type);	 
	  
-- create clob from remote
SQL> create table remote ( x varchar2(10), y clob );
Table created.

SQL> insert into remote values ( 'hello', 'This is a test');
1 row created.

SQL> commit;
Commit complete.

Then I create a Temporary Table.

SQL> create global temporary table gtt on commit delete rows as select y from remote@dblink;
Table created.

SQL> select * from gtt;
no rows selected

SQL> desc gtt
 Name                                      Null?    Type
 ----------------------------------------- -------- ----------------------------
 Y                                                  CLOB

I populate the temporary table (local) with the CLOB data from the remote database.

SQL> insert into gtt select y from remote@dblink;
1 row created.

Now, I have the CLOB data locally, at my temporary table.

SQL> select * from gtt;

Y
--------------------------------------------------------------------------------
This is a test

connect to sqlplus 
sqlplus "user/password@(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(HOST=uklpdudasa-scan.uk.standardchartered.com)(PORT=1622))(CONNECT_DATA=(SERVER=DEDICATED)(SERVICE_NAME=LC_TSAAS_APP_UAT.uk.standardchartered.com)))"

sqlplus "user/pass@(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(Host=hostname.network)(Port=1521))(CONNECT_DATA=(SID=remote_SID)))"

step 1 - go to C:\Oracle\product\11.2.0\client_1\BIN 
step 2 - run --> sqlplus "LC_CMR_EIM_DEV_01/LC_CMR_EIM_DEV_01_123@hklpdddasb-scan.hk.standardchartered.com:1622/LC_CMREIM_DEV.hk.standardchartered.com"
 

sqlplus "LC_CMR_EIM_DEV_01/LC_CMR_EIM_DEV_01_123@(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(Host=hklpdddasb-scan.hk.standardchartered.com)(Port=1622))(CONNECT_DATA=(SID=LC_CMREIM_DEV.hk.standardchartered.com)))"


sqlplus "LC_CMR_EIM_DEV_01/LC_CMR_EIM_DEV_01_123@(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(Host=hklpdddasb-scan.hk.standardchartered.com)(Port=1622))(CONNECT_DATA=(SID=LC_CMREIM_DEV)))"

sqlplus "LC_CMR_EIM_DEV_01/LC_CMR_EIM_DEV_01_123@hklpdddasb-scan.hk.standardchartered.com:1622/LC_CMREIM_DEV.hk.standardchartered.com"

-- run host command
CREATE OR REPLACE FUNCTION exec_host_command( lc_cmd IN VARCHAR2 )
RETURN INTEGER IS
ln_status NUMBER;
lc_errormsg VARCHAR2(80);
lc_pipe_name VARCHAR2(30);
BEGIN
lc_pipe_name := ‘HOST_PIPE’;
dbms_pipe.pack_message( lc_cmd );
ln_status := dbms_pipe.send_message(lc_pipe_name);
RETURN ln_status;
END;
/

postgresql
# psql -d cash_aml --> will connect to cash_aml database
cash_aml=# \dt *. --> List all tables
cash_aml=# \dv *. --> List all views
cash_aml=# \d table_name --> describe table name
cash_aml=# \d+ view_name --> describe and view definition (sql)
cash_aml=# \h delete --> display help on "delete"
cash_aml=# \i basic.sql --> executes all sql commands from basic.sql
cash_aml=# ALTER TABLE products RENAME COLUMN product_no TO product_number; -- rename column
cash_aml=# ALTER TABLE products RENAME TO items;  -- rename table
cash_aml=# ALTER TABLE products ALTER COLUMN price TYPE numeric(10,2); -- alter column datatype
cash_aml=# ALTER TABLE products ADD COLUMN description text; -- add column
cash_aml=# ALTER TABLE products DROP COLUMN description; -- drop column
cash_aml=# Copy (Select * From foo) To '/tmp/test.csv' With CSV; -- generate csv file -- only superuser
cash_aml=# Copy (Select * From foo) To '/tmp/test.txt' With TXT; -- generate txt file

AWK
1) print first column
awk -F "," '{print $1}' one.csv

Name	Description
A3	Accurate, Adaptable, and Accessible Error Metrics for Predictive Models
abbyyR	Access to Abbyy Optical Character Recognition (OCR) API
abc	Tools for Approximate Bayesian Computation (ABC)
ABCanalysis	Computed ABC Analysis
abc.data	Data Only: Tools for Approximate Bayesian Computation (ABC)
abcdeFBA	ABCDE_FBA: A-Biologist-Can-Do-Everything of Flux Balance Analysis with this package
ABCoptim	Implementation of Artificial Bee Colony (ABC) Optimization
ABCp2	Approximate Bayesian Computational Model for Estimating P2
abcrf	Approximate Bayesian Computation via Random Forests
abctools	Tools for ABC Analyses
abd	The Analysis of Biological Data
abf2	Load Gap-Free Axon ABF2 Files
ABHgenotypeR	Easy Visualization of ABH Genotypes
abind	Combine Multidimensional Arrays
abn	Modelling Multivariate Data with Additive Bayesian Networks
abodOutlier	Angle-Based Outlier Detection
AbsFilterGSEA	Improved False Positive Control of Gene-Permuting GSEA with Absolute Filtering
abundant	Abundant regression and high-dimensional principal fitted components
ACA	Abrupt Change-Point or Aberration Detection in Point Series
acc	Functions for Processing and Analyzing Accelerometer Data
accelerometry	Functions for Processing Minute-to-Minute Accelerometer Data
accelmissing	Missing Value Imputation for Accelerometer Data
AcceptanceSampling	Creation and Evaluation of Acceptance Sampling Plans
ACCLMA	ACC & LMA Graph Plotting
accrual	Bayesian Accrual Prediction
accrued	Data Quality Visualization Tools for Partially Accruing Data
ACD	Categorical data analysis with complete or missing responses
ACDm	Tools for Autoregressive Conditional Duration Models
acepack	ace() and avas() for selecting regression transformations
ACEt	Estimating Age Modification Effect on Genetic and Environmental Variance Components in Twin Models
acid	Analysing Conditional Income Distributions
acm4r	Align-and-Count Method comparisons of RFLP data
acmeR	Implements ACME Estimator of Bird and Bat Mortality by Wind Turbines
ACNE	Affymetrix SNP Probe-Summarization using Non-Negative Matrix Factorization
acnr	Annotated Copy-Number Regions
acopula	Modelling dependence with multivariate Archimax (or any user-defined continuous) copulas
acp	Autoregressive Conditional Poisson
aCRM	Convenience functions for analytical Customer Relationship Management
acs	Download, Manipulate, and Present American Community Survey and Decennial Data from the US Census
ACSNMineR	Gene Enrichment Analysis from ACSN Maps or GMT Files
acss	Algorithmic Complexity for Short Strings
acss.data	Data Only: Algorithmic Complexity of Short Strings (Computed via Coding Theorem Method)
ACSWR	A Companion Package for the Book "A Course in Statistics with R"
ACTCD	Asymptotic Classification Theory for Cognitive Diagnosis
Actigraphy	Actigraphy Data Analysis
activity	Animal Activity Statistics
activpalProcessing	Process activPAL Events Files
actuar	Actuarial Functions and Heavy Tailed Distributions
ActuDistns	Functions for actuarial scientists
ada	ada: an R package for stochastic boosting
adabag	Applies Multiclass AdaBoost.M1, SAMME and Bagging
adagio	Discrete and Global Optimization Routines
AdapEnetClass	A Class of Adaptive Elastic Net Methods for Censored Data
adaptDA	Adaptive Mixture Discriminant Analysis
AdaptFit	Adaptive Semiparametic Regression
AdaptFitOS	Adaptive Semiparametric Regression with Simultaneous Confidence Bands
AdaptGauss	Gaussian Mixture Models (GMM)
AdaptiveSparsity	Adaptive Sparsity Models
adaptivetau	Tau-leaping stochastic simulation
adaptMCMC	Implementation of a generic adaptive Monte Carlo Markov Chain sampler
adaptsmoFMRI	Adaptive Smoothing of FMRI Data
adaptTest	Adaptive two-stage tests
addhaz	Binomial and Multinomial Additive Hazards Models
addhazard	Fit Additive Hazards Models for Survival Analysis
additivityTests	Additivity Tests in the Two Way Anova with Single Sub-class Numbers
addreg	Additive Regression for Discrete Data
ADDT	A Package for Analysis of Accelerated Destructive Degradation Test Data
ade4	Analysis of Ecological Data : Exploratory and Euclidean Methods in Environmental Sciences
ade4TkGUI	'ade4' Tcl/Tk Graphical User Interface
adegenet	Exploratory Analysis of Genetic and Genomic Data
adegraphics	An S4 Lattice-Based Package for the Representation of Multivariate Data
adehabitat	Analysis of Habitat Selection by Animals
adehabitatHR	Home Range Estimation
adehabitatHS	Analysis of Habitat Selection by Animals
adehabitatLT	Analysis of Animal Movements
adehabitatMA	Tools to Deal with Raster Maps
adephylo	adephylo: exploratory analyses for the phylogenetic comparative method
AdequacyModel	Adequacy of probabilistic models and generation of pseudo-random numbers
ADGofTest	Anderson-Darling GoF test
adhoc	calculate ad hoc distance thresholds for DNA barcoding identification
adimpro	Adaptive Smoothing of Digital Images
adlift	An adaptive lifting scheme algorithm
ADM3	An Interpretation of the ADM method - automated detection algorithm
AdMit	Adaptive Mixture of Student-t distributions
ADMMnet	Regularized Model with Selecting the Number of Non-Zeros
ADPclust	Fast Clustering Using Adaptive Density Peak Detection
ads	Spatial point patterns analysis
adwave	Wavelet Analysis of Genomic Data from Admixed Populations
AEDForecasting	Change Point Analysis in ARIMA Forecasting
aemo	Download and process AEMO price and demand data
AER	Applied Econometrics with R
AF	Model-Based Estimation of Confounder-Adjusted Attributable Fractions
afex	Analysis of Factorial Experiments
AFLPsim	Hybrid Simulation and Genome Scan for Dominant Markers
AFM	Atomic Force Microscope Image Analysis
aftgee	Accelerated Failure Time Model with Generalized Estimating Equations
AGD	Analysis of Growth Data
AggregateR	Aggregate Numeric, Date and Categorical Variables by an ID
agop	Aggregation Operators and Preordered Sets
agRee	Various Methods for Measuring Agreement
Agreement	Statistical Tools for Measuring Agreement
agricolae	Statistical Procedures for Agricultural Research
agridat	Agricultural Datasets
agrmt	Calculate Agreement or Consensus in Ordered Rating Scales
AGSDest	Estimation in Adaptive Group Sequential Trials
agsemisc	Miscellaneous plotting and utility functions
ahaz	Regularization for semiparametric additive hazards regression
ahp	Analytic Hierarchy Process
AHR	Estimation and Testing of Average Hazard Ratios
AICcmodavg	Model Selection and Multimodel Inference Based on (Q)AIC(c)
AID	Estimation of Box-Cox Power Transformation Parameter
aidar	Tools for reading AIDA (http://aida.freehep.org/) files into R
AIM	AIM: adaptive index model
Ake	Associated Kernel Estimations
akima	Interpolation of Irregularly and Regularly Spaced Data
akmeans	Adaptive Kmeans algorithm based on threshold
alabama	Constrained Nonlinear Optimization
alakazam	Immunoglobulin Clonal Lineage and Diversity Analysis
ald	The Asymmetric Laplace Distribution
ALDqr	Quantile Regression Using Asymmetric Laplace Distribution
aLFQ	Estimating Absolute Protein Quantities from Label-Free LC-MS/MS Proteomics Data
AlgDesign	Algorithmic Experimental Design
AlgebraicHaploPackage	Haplotype Two Snips Out of a Paired Group of Patients
algstat	Algebraic statistics in R
AlignStat	Comparison of Alternative Multiple Sequence Alignments
alineR	Alignment of Phonetic Sequences Using the 'ALINE' Algorithm
ALKr	Generate Age-Length Keys for fish populations
allan	Automated Large Linear Analysis Node
allanvar	Allan Variance Analysis
alleHap	Allele Imputation and Haplotype Reconstruction from Pedigree Databases
allelematch	Identifying unique multilocus genotypes where genotyping error and missing data may be present
AlleleRetain	Allele Retention, Inbreeding, and Demography
allelic	A fast, unbiased and exact allelic exact test
AllPossibleSpellings	Computes all of a word's possible spellings
alm	R Client for the Lagotto Altmetrics Platform
alphahull	Generalization of the Convex Hull of a Sample of Points in the Plane
alphaOutlier	Obtain Alpha-Outlier Regions for Well-Known Probability Distributions
alphashape3d	Implementation of the 3D Alpha-Shape for the Reconstruction of 3D Sets from a Point Cloud
alr3	Data to accompany Applied Linear Regression 3rd edition
alr4	Data to accompany Applied Linear Regression 4rd edition
ALS	Multivariate Curve Resolution Alternating Least Squares (MCR-ALS)
ALSCPC	Accelerated line search algorithm for simultaneous orthogonal transformation of several positive definite symmetric matrices to nearly diagonal form
altmeta	Alternative Meta-Analysis Methods
ALTopt	Optimal Experimental Designs for Accelerated Life Testing
amap	Another Multidimensional Analysis Package
AMAP.Seq	Compare Gene Expressions from 2-Treatment RNA-Seq Experiments
AMCP	A Multiple Comparison Perspective
ameco	European Commission Annual Macro-Economic (AMECO) Database
amei	Adaptive Management of Epidemiological Interventions
Amelia	A Program for Missing Data
amen	Additive and Multiplicative Effects Models for Networks and Relational Data
AmericanCallOpt	This package includes pricing function for selected American call options with underlying assets that generate payouts
AMGET	Post-processing tool for ADAPT 5
aml	Adaptive Mixed LASSO
AMOEBA	A Multidirectional Optimum Ecotope-Based Algorithm
AMORE	A MORE flexible neural network package
AmpliconDuo	Statistical Analysis of Amplicon Data of the Same Sample to Identify Artefacts
anacor	Simple and Canonical Correspondence Analysis
analogsea	Interface to 'Digital Ocean'
analogue	Analogue and Weighted Averaging Methods for Palaeoecology
analogueExtra	Additional functions for use with the analogue package
analyz	Model Layer for Automatic Data Analysis via CSV File Interpretation
AnalyzeFMRI	Functions for analysis of fMRI datasets stored in the ANALYZE or NIFTI format
AnalyzeTS	Analyze Time Series
anametrix	Bidirectional connector to Anametrix API
anapuce	Tools for microarray data analysis
AncestryMapper	Assigning Ancestry Based on Population References
anchors	Statistical analysis of surveys with anchoring vignettes
AnDE	An extended Bayesian Learning Technique developed by Dr. Geoff Webb
andrews	Andrews curves
anesrake	ANES Raking Implementation
anfis	Adaptive Neuro Fuzzy Inference System in R
AnglerCreelSurveySimulation	Simulate a Bus Route Creel Survey of Anglers
animalTrack	Animal track reconstruction for high frequency 2-dimensional (2D) or 3-dimensional (3D) movement data
animation	A Gallery of Animations in Statistics and Utilities to Create Animations
anim.plots	Simple Animated Plots For R
AnnotLists	AnnotLists: A tool to annotate multiple lists from a specific annotation file
anoint	Analysis of Interactions
ANOM	Analysis of Means
anominate	alpha-NOMINATE Ideal Point Estimator
anonymizer	Anonymize Data Containing Personally Identifiable Information
AnthropMMD	A GUI for Mean Measures of Divergence
Anthropometry	Statistical Methods for Anthropometric Data
antitrust	Tools for Antitrust Practitioners
AntWeb	programmatic interface to the AntWeb
aod	Analysis of Overdispersed Data
aods3	Analysis of Overdispersed Data using S3 methods
aoos	Another Object Orientation System
aop	Adverse Outcome Pathway Analysis
aoristic	aoristic analysis with spatial output (kml)
ApacheLogProcessor	Process the Apache Web Server Log Combined Files
apaStyle	Generate APA Tables for MS Word
apaTables	Create American Psychological Association (APA) Style Tables
apc	Age-Period-Cohort Analysis
apcluster	Affinity Propagation Clustering
ape	Analyses of Phylogenetics and Evolution
apex	Phylogenetic Methods for Multiple Gene Data
aplore3	Datasets from Hosmer, Lemeshow and Sturdivant, "Applied Logistic Regression" (3rd ed.)
aplpack	Another Plot PACKage: stem.leaf, bagplot, faces, spin3R, plotsummary, plothulls, and some slider functions
apmsWAPP	Pre- and Postprocessing for AP-MS data analysis based on spectral counts
appell	Compute Appell's F1 hypergeometric function
apple	Approximate Path for Penalized Likelihood Estimators
AppliedPredictiveModeling	Functions and Data Sets for 'Applied Predictive Modeling'
appnn	Amyloid Propensity Prediction Neural Network
approximator	Bayesian prediction of complex computer codes
aprean3	Datasets from Draper and Smith "Applied Regression Analysis" (3rd Ed., 1998)
apricom	Tools for the a Priori Comparison of Regression Modelling Strategies
aprof	Amdahl's Profiler, Directed Optimization Made Easy
APSIM	General Utility Functions for the 'Agricultural Production Systems Simulator'
APSIMBatch	Analysis the output of Apsim software
apsimr	Edit, Run and Evaluate APSIM Simulations Easily Using R
apsrtable	apsrtable model-output formatter for social science
apt	Asymmetric Price Transmission
apTreeshape	Analyses of Phylogenetic Treeshape
aqfig	Functions to help display air quality model output and monitoring data
aqp	Algorithms for Quantitative Pedology
aqr	Interface methods to use with an ActiveQuant Master Server
AquaEnv	AquaEnv - an integrated development toolbox for aquatic chemical model generation
AR1seg	Segmentation of an autoregressive Gaussian process of order 1
ArArRedux	Rigorous Data Reduction and Error Propagation of Ar40 / Ar39 Data
archdata	Example Datasets from Archaeological Research
archetypes	Archetypal Analysis
archiDART	Plant Root System Architecture Analysis Using DART and RSML Files
archivist	Tools for Storing, Restoring and Searching for R Objects
archivist.github	Tools for Archiving, Managing and Sharing R Objects via GitHub
ArDec	Time series autoregressive-based decomposition
arf3DS4	Activated Region Fitting, fMRI data analysis (3D)
arfima	Fractional ARIMA (and Other Long Memory) Time Series Modeling
ArfimaMLM	Arfima-MLM Estimation For Repeated Cross-Sectional Data
argosfilter	Argos locations filter
argparse	Command line optional and positional argument parser
argparser	Command-Line Argument Parser
ArgumentCheck	Improved Communication to Users with Respect to Problems in Function Arguments
arm	Data Analysis Using Regression and Multilevel/Hierarchical Models
arnie	"Arnie" box office records 1982-2014
aroma.affymetrix	Analysis of Large Affymetrix Microarray Data Sets
aroma.apd	A Probe-Level Data File Format Used by 'aroma.affymetrix' [deprecated]
aroma.cn	Copy-Number Analysis of Large Microarray Data Sets
aroma.core	Core Methods and Classes Used by 'aroma.*' Packages Part of the Aroma Framework
ARPobservation	Tools for Simulating Direct Behavioral Observation Recording Procedures Based on Alternating Renewal Processes
aRpsDCA	Arps Decline Curve Analysis in R
arrApply	Apply a Function to a Margin of an Array
ArrayBin	Binarization of numeric data arrays
arrayhelpers	Convenience functions for arrays
ars	Adaptive Rejection Sampling
ART	Aligned Rank Transform for Nonparametric Factorial Analysis
artfima	ARTFIMA Model Estimation
ARTIVA	Time-Varying DBN Inference with the ARTIVA (Auto Regressive TIme VArying) Model
ARTool	Aligned Rank Transform
ARTP	Gene and Pathway p-values computed using the Adaptive Rank Truncated Product
ARTP2	Pathway and Gene-Level Association Test
arules	Mining Association Rules and Frequent Itemsets
arulesNBMiner	Mining NB-Frequent Itemsets and NB-Precise Rules
arulesSequences	Mining Frequent Sequences
arulesViz	Visualizing Association Rules and Frequent Itemsets
aRxiv	Interface to the arXiv API
asaur	Data Sets for "Applied Survival Analysis Using R""
asbio	A Collection of Statistical Tools for Biologists
ascii	Export R objects to several markup languages
asd	Simulations for adaptive seamless designs
asdreader	Reading ASD Binary Files in R
ash	David Scott's ASH Routines
asht	Applied Statistical Hypothesis Tests
AsioHeaders	'Asio' C++ Header Files
ASMap	Linkage Map Construction using the MSTmap Algorithm
asnipe	Animal Social Network Inference and Permutations for Ecologists
aspace	A collection of functions for estimating centrographic statistics and computational geometries for spatial point patterns
ASPBay	Bayesian Inference on Causal Genetic Variants using Affected Sib-Pairs Data
aspect	A General Framework for Multivariate Analysis with Optimal Scaling
aSPU	Adaptive Sum of Powered Score Test
asremlPlus	Augments the Use of 'Asreml' in Fitting Mixed Models
assertive	Readable Check Functions to Ensure Code Integrity
assertive.base	A Lightweight Core of the 'assertive' Package
assertive.code	Assertions to Check Properties of Code
assertive.data	Assertions to Check Properties of Data
assertive.data.uk	Assertions to Check Properties of Strings
assertive.data.us	Assertions to Check Properties of Strings
assertive.datetimes	Assertions to Check Properties of Dates and Times
assertive.files	Assertions to Check Properties of Files
assertive.matrices	Assertions to Check Properties of Matrices
assertive.models	Assertions to Check Properties of Models
assertive.numbers	Assertions to Check Properties of Numbers
assertive.properties	Assertions to Check Properties of Variables
assertive.reflection	Assertions for Checking the State of R
assertive.sets	Assertions to Check Properties of Sets
assertive.strings	Assertions to Check Properties of Strings
assertive.types	Assertions to Check Types of Variables
assertr	Assertive Programming for R Analysis Pipelines
assertthat	Easy pre and post assertions
AssetPricing	Optimal pricing of assets with fixed expiry date
assist	A Suite of R Functions Implementing Spline Smoothing Techniques
AssocTests	Genetic Association Studies
assortnet	Calculate the Assortativity Coefficient of Weighted and Binary Networks
AssotesteR	Statistical Tests for Genetic Association Studies
aster	Aster Models
aster2	Aster Models
astro	Astronomy Functions, Tools and Routines
astrochron	A Computational Tool for Astrochronology
astrodatR	Astronomical Data
astroFns	Astronomy: time and position functions, misc. utilities
astrolibR	Astronomy Users Library
astsa	Applied Statistical Time Series Analysis
asVPC	Average Shifted Visual Predictive Checks
asymLD	Asymmetric Linkage Disequilibrium (ALD) for Polymorphic Genetic Data
asympTest	Asymptotic statistic
AsynchLong	Regression Analysis of Sparse Asynchronous Longitudinal Data
asypow	Calculate Power Utilizing Asymptotic Likelihood Ratio Methods
ATE	Inference for Average Treatment Effects using Covariate Balancing
AtelieR	A GTK GUI for teaching basic concepts in statistical inference, and doing elementary bayesian tests
atmcmc	Automatically Tuned Markov Chain Monte Carlo
ATmet	Advanced Tools for Metrology
AtmRay	Acoustic Traveltime Calculations for 1-D Atmospheric Models
aTSA	Alternative Time Series Analysis
atsd	Support Querying Axibase Time-Series Database
attribrisk	Population Attributable Risk
AUC	Threshold independent performance measures for probabilistic classifiers
aucm	AUC Maximization
AUCRF	Variable Selection with Random Forest and the Area Under the Curve
audio	Audio Interface for R
audiolyzR	audiolyzR: Give your data a listen
audit	Bounds for Accounting Populations
auRoc	Various Methods to Estimate the AUC
autoencoder	Sparse Autoencoder for Automatic Learning of Representative Features from Unlabeled Data
automap	Automatic interpolation package
AutoModel	Automated Hierarchical Multiple Regression with Assumptions Checking
autopls	Partial Least Squares Regression with Backward Selection of Predictors
AutoregressionMDE	Minimum Distance Estimation in Autoregressive Model
AutoSEARCH	General-to-Specific (GETS) Modelling
autovarCore	Automated Vector Autoregression Models and Networks
averisk	Calculation of Average Population Attributable Fractions and Confidence Intervals
aws	Adaptive Weights Smoothing
awsMethods	Class and Methods Definitions for Packages 'aws', 'adimpro', 'fmri' and 'dti'
aws.signature	Amazon Web Services Request Signatures
aylmer	A generalization of Fisher's exact test
AzureML	Interface with Azure Machine Learning Datasets, Experiments and Web Services
B2Z	Bayesian Two-Zone Model
b6e6rl	Adaptive differential evolution, b6e6rl variant
babar	Bayesian Bacterial Growth Curve Analysis in R
babel	Ribosome profiling data analysis
BaBooN	Bayesian Bootstrap Predictive Mean Matching - Multiple and Single Imputation for Discrete Data
babynames	US Baby Names 1880-2014
BACA	Bubble Chart to Compare Biological Annotations by using DAVID
BacArena	Modeling Framework for Cellular Communities in their Environments
BACCO	Bayesian Analysis of Computer Code Output (BACCO)
backblazer	Bindings to the Backblaze B2 API
backpipe	Backward Pipe Operator
backports	Reimplementations of Functions Introduced Since R-3.0.0
backShift	Learning Causal Cyclic Graphs from Unknown Shift Interventions
backtest	Exploring Portfolio-Based Conjectures About Financial Instruments
backtestGraphics	Interactive Graphics for Portfolio Data
BACprior	Choice of the Hyperparameter Omega in the Bayesian Adjustment for Confounding (BAC) Algorithm
bacr	Bayesian Adjustment for Confounding
BAEssd	Bayesian Average Error approach to Sample Size Determination
Bagidis	BAses GIving DIStances
bagRboostR	Ensemble bagging and boosting classifiers
BalancedSampling	Balanced and Spatially Balanced Sampling
BaM	Functions and datasets for books by Jeff Gill
bamdit	Bayesian Meta-Analysis of Diagnostic Test Data
BAMMtools	Analysis and Visualization of Macroevolutionary Dynamics on Phylogenetic Trees
bandit	Functions for simple A/B split test and multi-armed bandit analysis
BANFF	Bayesian Network Feature Finder
BANOVA	Hierarchical Bayesian ANOVA Models
bapred	Batch Effect Removal (in Phenotype Prediction using Gene Data)
barcode	Barcode distribution plots
Barnard	Barnard's Unconditional Test
bartMachine	Bayesian Additive Regression Trees
bartMachineJARs	bartMachine JARs
BAS	Bayesian Model Averaging using Bayesian Adaptive Sampling
base64	Base 64 encoder/decoder
base64enc	Tools for base64 encoding
basefun	Infrastructure for Computing with Basis Functions
baseline	Baseline Correction of Spectra
basicspace	Recovering a Basic Space from Issue Scales
BASIX	BASIX: An efficient C/C++ toolset for R
BaSTA	Age-Specific Survival Analysis from Incomplete Capture-Recapture/Recovery Data
BAT	Biodiversity Assessment Tools
batade	HTML reports and so on
batch	Batching Routines in Parallel and Passing Command-Line Arguments to R
BatchExperiments	Statistical Experiments on Batch Computing Clusters
BatchJobs	Batch Computing with R
batchmeans	Consistent Batch Means Estimation of Monte Carlo Standard Errors
batman	Convert Categorical Representations of Logicals to Actual Logicals
batteryreduction	An R Package for Data Reduction by Battery Reduction
BayClone2	Bayesian Feature Allocation Model for Tumor Heterogeneity
BayesBD	Bayesian Boundary Detection in Images
bayesboot	An Implementation of Rubin's (1981) Bayesian Bootstrap
BayesBridge	Bridge Regression
BayesComm	Bayesian Community Ecology Analysis
bayescount	Power Calculations and Bayesian Analysis of Count Distributions and FECRT Data using MCMC
BayesCR	Bayesian Analysis of Censored Regression Models Under Scale Mixture of Skew Normal Distributions
BayesDA	Functions and Datasets for the book "Bayesian Data Analysis"
bayesDccGarch	The Bayesian Dynamic Conditional Correlation GARCH Model
bayesDem	Graphical User Interface for bayesTFR, bayesLife and bayesPop
BayesFactor	Computation of Bayes Factors for Common Designs
bayesGARCH	Bayesian Estimation of the GARCH(1,1) Model with Student-t Innovations
bayesGDS	Scalable Rejection Sampling for Bayesian Hierarchical Models
BayesGESM	Bayesian Analysis of Generalized Elliptical Semi-Parametric Models and Flexible Measurement Error Models
BayesianAnimalTracker	Bayesian Melding of GPS and DR Path for Animal Tracking
Bayesianbetareg	Bayesian Beta regression: joint mean and precision modeling
BayesLCA	Bayesian Latent Class Analysis
bayesLife	Bayesian Projection of Life Expectancy
BayesLogit	Logistic Regression
bayesm	Bayesian Inference for Marketing/Micro-Econometrics
BayesMAMS	Designing Bayesian Multi-Arm Multi-Stage Studies
bayesMCClust	Mixtures-of-Experts Markov Chain Clustering and Dirichlet Multinomial Clustering
BayesMed	Default Bayesian Hypothesis Tests for Correlation, Partial Correlation, and Mediation
bayesmeta	Bayesian Random-Effects Meta-Analysis
bayesmix	Bayesian Mixture Models with JAGS
BayesMixSurv	Bayesian Mixture Survival Models using Additive Mixture-of-Weibull Hazards, with Lasso Shrinkage and Stratification
BayesNI	BayesNI: Bayesian Testing Procedure for Noninferiority with Binary Endpoints
bayesPop	Probabilistic Population Projection
bayespref	Hierarchical Bayesian analysis of ecological count data
bayesQR	Bayesian quantile regression
bayess	Bayesian Essentials with R
BayesSAE	Bayesian Analysis of Small Area Estimation
BayesSingleSub	Computation of Bayes factors for interrupted time-series designs
BayesSummaryStatLM	MCMC Sampling of Bayesian Linear Models via Summary Statistics
bayesSurv	Bayesian Survival Regression with Flexible Error and Random Effects Distributions
bayesTFR	Bayesian Fertility Projection
Bayesthresh	Bayesian thresholds mixed-effects models for categorical data
BayesTree	Bayesian Additive Regression Trees
BayesValidate	BayesValidate Package
BayesVarSel	Bayes Factors, Model Choice and Variable Selection in Linear Models
BayesX	R Utilities Accompanying the Software Package BayesX
BayesXsrc	R Package Distribution of the BayesX C++ Sources
BayHap	Bayesian analysis of haplotype association using Markov Chain Monte Carlo
BayHaz	R Functions for Bayesian Hazard Rate Estimation
BaylorEdPsych	R Package for Baylor University Educational Psychology Quantitative Courses
bayou	Bayesian Fitting of Ornstein-Uhlenbeck Models to Phylogenies
BaySIC	Bayesian Analysis of Significantly Mutated Genes in Cancer
BAYSTAR	On Bayesian analysis of Threshold autoregressive model (BAYSTAR)
BB	Solving and Optimizing Large-Scale Nonlinear Systems
bbefkr	Bayesian bandwidth estimation and semi-metric selection for the functional kernel regression with unknown error density
bbemkr	Bayesian bandwidth estimation for multivariate kernel regression with Gaussian error
BBEST	Bayesian Estimation of Incoherent Neutron Scattering Backgrounds
BBmisc	Miscellaneous Helper Functions for B. Bischl
bbmle	Tools for General Maximum Likelihood Estimation
BBMM	Brownian bridge movement model
bbo	Biogeography-Based Optimization
BBRecapture	Bayesian Behavioural Capture-Recapture Models
bc3net	Gene Regulatory Network Inference with Bc3net
BCA	Business and Customer Analytics
BCBCSF	Bias-Corrected Bayesian Classification with Selected Features
BCDating	Business Cycle Dating and Plotting Tools
BcDiag	Diagnostics Plots for Bicluster Data
BCE	Bayesian composition estimator: estimating sample (taxonomic) composition from biomarker data
BCEA	Bayesian Cost Effectiveness Analysis
BCEE	The Bayesian Causal Effect Estimation Algorithm
BCEs0	Bayesian Models for Cost-Effectiveness Analysis in the Presence of Structural Zero Costs
Bchron	Radiocarbon Dating, Age-Depth Modelling, Relative Sea Level Rate Estimation, and Non-Parametric Phase Modelling
Bclim	Bayesian Palaeoclimate Reconstruction from Pollen
bclust	Bayesian Hierarchical Clustering Using Spike and Slab Models
bcp	Bayesian Analysis of Change Point Problems
bcpa	Behavioral change point analysis of animal movement
bcpmeta	Bayesian Multiple Changepoint Detection Using Metadata
BCRA	Breast Cancer Risk Assessment
bcRep	Advanced Analysis of B Cell Receptor Repertoire Data
bcrm	Bayesian Continual Reassessment Method for Phase I Dose-Escalation Trials
bcrypt	'Blowfish' Password Hashing Algorithm
bcv	Cross-Validation for the SVD (Bi-Cross-Validation)
bda	Density Estimation for Grouped Data
bde	Bounded Density Estimation
BDgraph	Bayesian Graph Selection Based on Birth-Death MCMC Approach
bdots	Bootstrapped Differences of Time Series
bdpv	Inference and design for predictive values in binary diagnostic tests
bdrift	Beta Drift Analysis
bdscale	Remove Weekends and Holidays from ggplot2 Axes
bdsmatrix	Routines for Block Diagonal Symmetric matrices
bdvis	Biodiversity Data Visualizations
BDWreg	Bayesian Inference for Discrete Weibull Regression
bdynsys	Bayesian Dynamical System Model
beadarrayFilter	Bead filtering for Illumina bead arrays
beadarrayMSV	Analysis of Illumina BeadArray SNP data including MSV markers
beanplot	Visualization via Beanplots (like Boxplot/Stripchart/Violin Plot)
BEANSP	Bayesian Estimate of Age-specific Nest Survival Probabilities
BEDASSLE	Quantifies effects of geo/eco distance on genetic differentiation
BEDMatrix	Matrices Backed by Binary PED Files (PLINK)
bedr	Genomic Region Processing using Tools Such as Bedtools, Bedops and Tabix
beepr	Easily Play Notification Sounds on any Platform
beeswarm	The Bee Swarm Plot, an Alternative to Stripchart
benchden	28 benchmark densities from Berlinet/Devroye (1994)
benchmark	Benchmark Experiments Toolbox
Benchmarking	Benchmark and Frontier Analysis Using DEA and SFA
benchmarkme	Crowd Sourced System Benchmarks
benchmarkmeData	Data Set for the 'benchmarkme' Package
benford.analysis	Benford Analysis for Data Validation and Forensic Analytics
BenfordTests	Statistical Tests for Evaluating Conformity to Benford's Law
bentcableAR	Bent-Cable Regression for Independent Data or Autoregressive Time Series
BEQI2	Benthic Ecosystem Quality Index 2
ber	Batch Effects Removal
Bergm	Bayesian analysis for exponential random graph models
BerlinData	Easy access to Berlin related data
berryFunctions	Function Collection Related to Plotting and Hydrology
Bessel	Bessel – Bessel Functions Computations and Approximations
BEST	Bayesian Estimation Supersedes the t-Test
bestglm	Best Subset GLM
betafam	Detecting rare variants for quantitative traits using nuclear families
betalink	Beta-Diversity of Species Interactions
betapart	Partitioning beta diversity into turnover and nestedness components
betaper	Functions to incorporate taxonomic uncertainty on multivariate analyses of ecological data
betareg	Beta Regression
betas	Standardized Beta Coefficients
betategarch	Simulation, estimation and forecasting of Beta-Skew-t-EGARCH models
bethel	Bethel's algorithm
bezier	Bezier Curve and Spline Toolkit
bfa	Bayesian Factor Analysis
bfast	Breaks For Additive Season and Trend (BFAST)
bfork	Basic Unix Process Control
bfp	Bayesian Fractional Polynomials
bgeva	Binary Generalized Extreme Value Additive Models
bglm	Bayesian Estimation in Generalized Linear Models
BGLR	Bayesian Generalized Linear Regression
bgmm	Gaussian Mixture Modeling Algorithms And The Belief-based Mixture Modeling
BGPhazard	Markov Beta and Gamma Processes for Modeling Hazard Rates
BGSIMD	Block Gibbs Sampler with Incomplete Multinomial Distribution
BH	Boost C++ Header Files
Bhat	General likelihood exploration
BHH2	Useful Functions for Box, Hunter and Hunter II
biasbetareg	Bias correction of the parameter estimates of the beta regression model
BiasedUrn	Biased Urn Model Distributions
bibtex	bibtex parser
biclust	BiCluster Algorithms
BiDimRegression	Calculates the bidimensional regression between two 2D configurations
bifactorial	Inferences for bi- and trifactorial trial designs
BIFIEsurvey	Tools for Survey Statistics in Educational Assessment
bigalgebra	BLAS routines for native R matrices and big.matrix objects
biganalytics	Utilities for 'big.matrix' Objects from Package 'bigmemory'
bigdata	Big Data Analytics
BIGDAWG	Case-Control Analysis of Multi-Allelic Loci
bigGP	Distributed Gaussian Process Calculations
biglars	Scalable Least-Angle Regression and Lasso
biglasso	Big Lasso: Extending Lasso Model Fitting to Big Data in R
biglm	bounded memory linear and generalized linear models
bigmemory	Manage Massive Matrices with Shared Memory and Memory-Mapped Files
bigmemory.sri	A shared resource interface for Bigmemory Project packages
bigml	Bindings for the BigML API
bigpca	PCA, Transpose and Multicore Functionality for 'big.matrix' Objects
bigrquery	An Interface to Google's 'BigQuery' 'API'
bigRR	Generalized Ridge Regression (with special advantage for p >> n cases)
bigsplines	Smoothing Splines for Large Samples
bigtabulate	Table, Apply, and Split Functionality for Matrix and 'big.matrix' Objects
BigTSP	Top Scoring Pair based methods for classification
bild	BInary Longitudinal Data
bimetallic	Power for SNP analyses using silver standard cases
bimixt	Estimates Mixture Models for Case-Control Data
Binarize	Binarization of One-Dimensional Data
BinaryEMVS	Variable Selection for Binary Data Using the EM Algorithm
BinaryEPPM	Mean and Variance Modeling of Binary Data
binaryLogic	Binary Logic
binda	Multi-Class Discriminant Analysis using Binary Predictors
bindata	Generation of Artificial Binary Data
binequality	Methods for Analyzing Binned Income Data
bingat	Binary Graph Analysis Tools
binGroup	Evaluation and experimental design for binomial group testing
binhf	Haar-Fisz functions for binomial data
binMto	Asymptotic simultaneous confidence intervals for many-to-one comparisons of proportions
BinNonNor	Data Generation with Binary and Continuous Non-Normal Components
BinNor	Simultaneous generation of multivariate binary and normal variates
binom	Binomial Confidence Intervals For Several Parameterizations
binomen	'Taxonomic' Specification and Parsing Methods
binomialcftp	Generates binomial random numbers via the coupling from the past algorithm
binomlogit	Efficient MCMC for Binomial Logit Models
binomSamSize	Confidence intervals and sample size determination for a binomial proportion under simple random sampling and pooled sampling
binomTools	Performing diagnostics on binomial regression models
BinOrdNonNor	Concurrent Generation of Binary, Ordinal and Continuous Data
binr	Cut Numeric Values into Evenly Distributed Groups
binseqtest	Exact Binary Sequential Designs and Analysis
bio3d	Biological Structure Analysis
Biocomb	Feature Selection and Classification with the Embedded Validation Procedures for Biomedical Data Analysis
Biodem	Biodemography Functions
BiodiversityR	Package for Community Ecology and Suitability Analysis
BIOdry	Multilevel Modeling of Dendroclimatical Fluctuations
BioFTF	Biodiversity Assessment Using Functional Tools
biogas	Process Biogas Data and Predict Biogas Production
BioGeoBEARS	BioGeography with Bayesian (and Likelihood) Evolutionary Analysis in R Scripts
biogram	N-Gram Analysis of Biological Sequences
bioinactivation	Simulation of Dynamic Microbial Inactivation
bio.infer	Predict environmental conditions from biological observations
biom	An interface package (beta) for the BIOM file format
BioMark	Find Biomarkers in Two-Class Discrimination Problems
biomartr	Functional Annotation and Biological Data Retrieval with R
biomod2	Ensemble Platform for Species Distribution Modeling
BIOM.utils	Utilities for the BIOM (Biological Observation Matrix) Format
bionetdata	Biological and chemical data networks
BioPhysConnectoR	BioPhysConnectoR
bioPN	Simulation of deterministic and stochastic biochemical reaction networks using Petri Nets
biorxivr	Search and Download Papers from the bioRxiv Preprint Server
bios2mds	From BIOlogical Sequences to MultiDimensional Scaling
biosignalEMG	Tools for Electromyogram Signals (EMG) Analysis
BioStatR	Initiation à la Statistique avec R
biotools	Tools for Biometry and Applied Statistics in Agricultural Science
bipartite	Visualising bipartite networks and calculating some (ecological) indices
biplotbootGUI	Bootstrap on Classical Biplots and Clustering Disjoint Biplot
BiplotGUI	Interactive Biplots in R
BIPOD	BIPOD (Bayesian Inference for Partially Observed diffusions)
birdring	Methods to Analyse Ring Re-Encounter Data
birk	MA Birk's Functions
bisectr	Tools to find bad commits with git bisect
BiSEp	Toolkit to Identify Candidate Synthetic Lethality
bisoreg	Bayesian Isotonic Regression with Bernstein Polynomials
bit	A class for vectors of 1-bit booleans
bit64	A S3 Class for Vectors of 64bit Integers
bitops	Bitwise Operations
BiTrinA	Binarization and Trinarization of One-Dimensional Data
BivarP	Estimating the Parameters of Some Bivariate Distributions
bivarRIpower	Sample size calculations for bivariate longitudinal data
biwavelet	Conduct Univariate and Bivariate Wavelet Analyses
biwt	Functions to compute the biweight mean vector and covariance & correlation matrices
bizdays	Business Days Calculations and Utilities
BKPC	Bayesian Kernel Projection Classifier
BlakerCI	Blaker's Binomial Confidence Limits
BlandAltmanLeh	Plots (Slightly Extended) Bland-Altman Plots
blatr	Send Emails Using 'Blat' for Windows
Blaunet	Calculate and Analyze Blau Status for Measuring Social Distance
blavaan	Bayesian Latent Variable Analysis
BLCOP	Black-Litterman and Copula Opinion Pooling Frameworks
blender	Analyze biotic homogenization of landscapes
blighty	United Kingdom coastlines
blkergm	Fitting block ERGM given the block structure on social networks
blm	Binomial linear and linear-expit regression
blme	Bayesian Linear Mixed-Effects Models
blmeco	Data Files and Functions Accompanying the Book "Bayesian Data Analysis in Ecology using R, BUGS and Stan"
blockcluster	Coclustering Package for Binary, Categorical, Contingency and Continuous Data-Sets
blockmatrix	blockmatrix: Tools to solve algebraic systems with partitioned matrices
BlockMessage	Creates strings that show a text message in 8 by 8 block letters
blockmodeling	An R package for Generalized and classical blockmodeling of valued networks
blockmodels	Latent and Stochastic Block Model Estimation by a 'V-EM' Algorithm
blockrand	Randomization for block random clinical trials
blocksdesign	Nested Block Designs for Unstructured Treatments
blockseg	Two Dimensional Change-Points Detection
blockTools	Block, Assign, and Diagnose Potential Interference in Randomized Experiments
Blossom	Statistical Comparisons with Distance-Function Based Permutation Tests
BLR	Bayesian Linear Regression
blsAPI	Request Data from the U.S. Bureau of Labor Statistics API
BMA	Bayesian Model Averaging
bmd	Benchmark dose analysis for dose-response data
bmem	Mediation analysis with missing data using bootstrap
bmeta	Bayesian Meta-Analysis and Meta-Regression
BMhyd	PCM for Hybridization
Bmix	Bayesian Sampling for Stick-Breaking Mixtures
bmk	MCMC diagnostics package
bmmix	Bayesian multinomial mixture
BMN	The pseudo-likelihood method for pairwise binary markov networks
bmp	Read Windows Bitmap (BMP) images
bmrm	Bundle Methods for Regularized Risk Minimization Package
BMRV	Bayesian Models for Rare Variant Association Analysis
BMS	Bayesian Model Averaging Library
bnclassify	Learning Discrete Bayesian Network Classifiers from Data
BNDataGenerator	Data Generator based on Bayesian Network Model
bnlearn	Bayesian Network Structure Learning, Parameter Learning and Inference
bnnSurvival	Bagged k-Nearest Neighbors Survival Prediction
bnormnlr	Bayesian Estimation for Normal Heteroscedastic Nonlinear Regression Models
BNPdensity	Ferguson-Klass Type Algorithm for Posterior Normalized Random Measures
bnpmr	Bayesian monotonic nonparametric regression
BNPTSclust	A Bayesian Nonparametric Algorithm for Time Series Clustering
BNSP	Bayesian Non- and Semi-Parametric Model Fitting
bnstruct	Bayesian Network Structure Learning from Data with Missing Values
boa	Bayesian Output Analysis Program (BOA) for MCMC
bodenmiller	Profilling of Peripheral Blood Mononuclear Cells using CyTOF
BOG	Bacterium and Virus Analysis of Orthologous Groups (BOG) is a Package for Identifying Differentially Regulated Genes in the Light of Gene Functions
boilerpipeR	Interface to the Boilerpipe Java Library
BOIN	Bayesian Optimal INterval (BOIN) Design for Single-Agent and Drug- Combination Phase I Clinical Trials
bold	Interface to Bold Systems 'API'
Bolstad	Functions for Elementary Bayesian Inference
Bolstad2	Bolstad functions
BonEV	An Improved Multiple Testing Procedure for Controlling False Discovery Rates
boolean3	Boolean Binary Response Models
BoolNet	Construction, Simulation and Analysis of Boolean Networks
Boom	Bayesian Object Oriented Modeling
BoomSpikeSlab	MCMC for Spike and Slab Regression
boostmtree	Boosted Multivariate Trees for Longitudinal Data
boostr	A modular framework to bag or boost any estimation procedure
boostSeq	Optimized GWAS cohort subset selection for resequencing studies
boot	Bootstrap Functions (Originally by Angelo Canty for S)
bootES	Bootstrap Effect Sizes
bootLR	Bootstrapped Confidence Intervals for (Negative) Likelihood Ratio Tests
bootnet	Bootstrap Methods for Various Network Estimation Routines
BootPR	Bootstrap Prediction Intervals and Bias-Corrected Forecasting
bootRes	Bootstrapped Response and Correlation Functions
bootruin	A bootstrap test for the probability of ruin in the classical risk process
bootspecdens	Testing equality of spectral densities
bootsPLS	Bootstrap Subsamplings of Sparse Partial Least Squares - Discriminant Analysis for Classification and Signature Identification
bootStepAIC	Bootstrap stepAIC
bootstrap	Functions for the Book "An Introduction to the Bootstrap"
bootSVD	Fast, Exact Bootstrap Principal Component Analysis for High Dimensional Data
boottol	Bootstrap Tolerance Levels for Credit Scoring Validation Statistics
boral	Bayesian Ordination and Regression AnaLysis
Boruta	Wrapper Algorithm for All-Relevant Feature Selection
BoSSA	a Bunch of Structure and Sequence Analysis
boussinesq	Analytic Solutions for (ground-water) Boussinesq Equation
boxplotdbl	Double Box Plot for Two-Axes Correlation
boxr	Interface for the 'Box.com API'
bpca	Biplot of Multivariate Data Based on Principal Components Analysis
bpcp	Beta Product Confidence Procedure for Right Censored Data
bPeaks	bPeaks: an intuitive peak-calling strategy to detect transcription factor binding sites from ChIP-seq data in small eukaryotic genomes
bpkde	Back-Projected Kernel Density Estimation
bqtl	Bayesian QTL Mapping Toolkit
BradleyTerry2	Bradley-Terry Models
braidReports	Visualize Combined Action Response Surfaces and Report BRAID Analyses
braidrm	Fitting Dose Response with the BRAID Combined Action Model
BrailleR	Improved Access for Blind Users
brainGraph	Graph Theory Analysis of Brain MRI Data
brainR	Helper functions to misc3d and rgl packages for brain imaging
brainwaver	Basic wavelet analysis of multivariate time series with a visualisation and parametrisation using graph theory
braQCA	Bootstrapped Robustness Assessment for Qualitative Comparative Analysis
breakage	SICM pipette tip geometry estimation
breakaway	Species Richness Estimation and Modeling
breakpoint	An R Package for Multiple Break-Point Detection via the Cross-Entropy Method
bReeze	Functions for wind resource assessment
brew	Templating Framework for Report Generation
brewdata	Extracting Usable Data from the Grad Cafe Results Search
brglm	Bias reduction in binomial-response generalized linear models
bride	Brier score decomposition of probabilistic forecasts for binary events
brms	Bayesian Regression Models using Stan
brnn	Bayesian Regularization for Feed-Forward Neural Networks
Brobdingnag	Very large numbers in R
broman	Karl Broman's R Code
broom	Convert Statistical Analysis Objects into Tidy Data Frames
brotli	A Compression Format Optimized for the Web
brr	Bayesian Inference on the Ratio of Two Poisson Rates
brranching	Fetch 'Phylogenies' from Many Sources
BRugs	Interface to the 'OpenBUGS' MCMC Software
BSagri	Statistical methods for safety assessment in agricultural field trials
BSDA	Basic Statistics and Data Analysis
BSGS	Bayesian Sparse Group Selection
BSGW	Bayesian Survival Model with Lasso Shrinkage Using Generalized Weibull Regression
bshazard	Nonparametric Smoothing of the Hazard Function
BsMD	Bayes Screening and Model Discrimination
bspec	Bayesian Spectral Inference
bspmma	bspmma: Bayesian Semiparametric Models for Meta-Analysis
BSquare	Bayesian Simultaneous Quantile Regression
BSSasymp	Asymptotic Covariance Matrices of Some BSS Mixing and Unmixing Matrix Estimates
bssn	Birnbaum-Saunders Model Based on Skew-Normal Distribution
bst	Gradient Boosting
bsts	Bayesian Structural Time Series
btergm	Temporal Exponential Random Graph Models by Bootstrapped Pseudolikelihood
btf	Estimates univariate function via Bayesian trend filtering
BTLLasso	Modelling Heterogeneity in Paired Comparison Data
BTSPAS	Bayesian Time-Strat. Population Analysis
BTYD	Implementing Buy 'Til You Die Models
bujar	Buckley-James Regression for Survival Data with High-Dimensional Covariates
BurStFin	Burns Statistics Financial
BurStMisc	Burns Statistics miscellaneous
bursts	Markov model for bursty behavior in streams
bvarsv	Bayesian Analysis of a Vector Autoregressive Model with Stochastic Volatility and Time-Varying Parameters
bvenn	A Simple alternative to proportional Venn diagrams
bvls	The Stark-Parker algorithm for bounded-variable least squares
bvpSolve	Solvers for Boundary Value Problems of Differential Equations
BVS	Bayesian Variant Selection: Bayesian Model Uncertainty Techniques for Genetic Association Studies
bWGR	Bagging Whole-Genome Regression
c060	Extended Inference for Lasso and Elastic-Net Regularized Cox and Generalized Linear Models
c3net	Infering large-scale gene networks with C3NET
C50	C5.0 Decision Trees and Rule-Based Models
ca	Simple, Multiple and Joint Correspondence Analysis
cabootcrs	Bootstrap Confidence Regions for Correspondence Analysis
cacIRT	Classification Accuracy and Consistency under Item Response Theory
CaDENCE	Conditional Density Estimation Network Construction and Evaluation
CADFtest	This package performs the CADF unit root test proposed in Hansen (1995)
cAIC4	Conditional Akaike information criterion for lme4
Cairo	R graphics device using cairo graphics library for creating high-quality bitmap (PNG, JPEG, TIFF), vector (PDF, SVG, PostScript) and display (X11 and Win32) output
cairoDevice	Embeddable Cairo Graphics Device Driver
calACS	Calculations for All Common Subsequences
CALF	Coarse Approximation Linear Function
CALIBERrfimpute	Multiple imputation using MICE and Random Forest
calibrar	Automated Parameter Estimation for Complex (Ecological) Models
calibrate	Calibration of Scatterplot and Biplot Axes
calibrator	Bayesian calibration of complex computer codes
calmate	Improved Allele-Specific Copy Number of SNP Microarrays for Downstream Segmentation
CAM	Causal Additive Model (CAM)
CAMAN	Finite Mixture Models and Meta-Analysis Tools - Based on C.A.MAN
camel	Calibrated Machine Learning
camtrapR	Camera Trap Data Management and Preparation of Occupancy and Spatial Capture-Recapture Analyses
cancerTiming	Estimation of Temporal Ordering of Cancer Abnormalities
candisc	Visualizing Generalized Canonical Discriminant and Canonical Correlation Analysis
Canopy	Accessing Intra-Tumor Heterogeneity and Tracking Longitudinal and Spatial Clonal Evolutionary History by Next-Generation Sequencing
CANSIM2R	Directly Extracts Complete CANSIM Data Tables
cape	Combined analysis of pleiotropy and epistasis
caper	Comparative Analyses of Phylogenetics and Evolution in R
capm	Companion Animal Population Management
captioner	Numbers Figures and Creates Simple Captions
captr	Client for the Captricity API
capushe	CAlibrating Penalities Using Slope HEuristics
capwire	Estimates population size from non-invasive sampling
car	Companion to Applied Regression
CARBayes	Spatial Generalised Linear Mixed Models for Areal Unit Data
CARBayesdata	Data Sets Used in the Vignette Accompanying the CARBayes Package
CARBayesST	Spatio-Temporal Generalised Linear Mixed Models for Areal Unit Data
carcass	Estimation of the Number of Fatalities from Carcass Searches
cardidates	Identification of Cardinal Dates in Ecological Time Series
cardioModel	Cardiovascular Safety Exposure-Response Modeling in Early-Phase Clinical Studies
care	High-Dimensional Regression and CAR Score Variable Selection
CARE1	Statistical package for population size estimation in capture-recapture models
caret	Classification and Regression Training
caretEnsemble	Ensembles of Caret Models
caribou	Estimation of caribou abundance based on large scale aggregations monitored by radio telemetry
CarletonStats	Functions For Statistics Classes At Carleton College
CARLIT	Ecological Quality Ratios Calculation and Plot
caroline	A Collection of Database, Data Structure, Visualization, and Utility Functions for R
caRpools	CRISPR AnalyzeR for Pooled CRISPR Screens
CARrampsOcl	Reparameterized and marginalized posterior sampling for conditional autoregressive models, OpenCL implementation
cartography	Thematic Cartography
carx	Censored Autoregressive Model with Exogenous Covariates
caschrono	Séries temporelles avec R - Méthodes et cas
caseMatch	Identify Similar Cases for Qualitative Case Studies
cat	Analysis of categorical-variable datasets with missing values
catdap	Categorical Data Analysis Program Package
catdata	Categorical Data
CatDyn	Fishery Stock Assessment by Generalized Depletion Models
cate	High Dimensional Factor Analysis and Confounder Adjusted Testing and Estimation
catenary	Fits a Catenary to Given Points
CateSelection	Categorical Variable Selection Methods
cati	Community Assembly by Traits: Individuals and Beyond
catIrt	An R Package for Simulating IRT-Based Computerized Adaptive Tests
catnet	Categorical Bayesian Network Inference
caTools	Tools: moving window statistics, GIF, Base64, ROC AUC, etc
catR	Procedures to Generate Patterns under Computerized Adaptive Testing
catspec	Special models for categorical variables
causaldrf	Tools for Estimating Causal Dose Response Functions
causaleffect	Deriving Expressions of Joint Interventional Distributions and Transport Formulas in Causal Models
CausalFX	Methods for Estimating Causal Effects from Observational Data
CausalGAM	Estimation of Causal Effects with Generalized Additive Models
causalsens	Selection Bias Approach to Sensitivity Analysis for Causal Effects
Causata	Analysis utilities for binary classification and Causata users
CAvariants	Correspondence Analysis Variants
cba	Clustering for Business Analytics
CBPS	Covariate Balancing Propensity Score
cbsodataR	Statistics Netherlands (CBS) Open Data API Client
CCA	Canonical correlation analysis
CCAGFA	Bayesian Canonical Correlation Analysis and Group Factor Analysis
ccaPP	(Robust) Canonical Correlation Analysis via Projection Pursuit
cccd	Class Cover Catch Digraphs
ccChooser	Developing a core collections
cccp	Cone Constrained Convex Problems
cccrm	Concordance Correlation Coefficient for Repeated (and Non-Repeated) Measures
ccda	Combined Cluster and Discriminant Analysis
ccgarch	Conditional Correlation GARCH models
cchs	Cox Model for Case-Cohort Data with Stratified Subcohort-Selection
cclust	Convex Clustering Methods and Clustering Indexes
CCM	Correlation classification method (CCM)
CCMnet	Simulate Congruence Class Model for Networks
CCP	Significance Tests for Canonical Correlation Analysis (CCA)
CCpop	One and two locus GWAS of binary phenotype with case-control-population design
CCTpack	Cultural Consensus Theory applications to data
cda	Coupled dipole approximation of light scattering by clusters of nanoparticles
cdb	Reading and Writing Constant DataBases
cdcfluview	Retrieve U.S. Flu Season Data from the CDC FluView Portal
cdcsis	Conditional Distance Correlation and Its Related Feature Screening Method
cdfquantreg	Quantile Regression for Random Variables on the Unit Interval
CDFt	Statistical downscaling through CDF-transform
CDLasso	Coordinate Descent Algorithms for Lasso Penalized L1, L2, and Logistic Regression
CDM	Cognitive Diagnosis Modeling
CDNmoney	Components of Canadian Monetary and Credit Aggregates
cdom	R Functions to Model CDOM Spectra
CDROM	Phylogenetically Classifies Retention Mechanisms of Duplicate Genes from Gene Expression Data
cds	Constrained Dual Scaling for Detecting Response Styles
CDVine	Statistical Inference of C- And D-Vine Copulas
CEC	Cross-Entropy Clustering
cec2005benchmark	Benchmark for the CEC 2005 Special Session on Real-Parameter Optimization
cec2013	Benchmark functions for the Special Session and Competition on Real-Parameter Single Objective Optimization at CEC-2013
CEGO	Combinatorial Efficient Global Optimization
celestial	Collection of Common Astronomical Conversion Routines and Functions
cellranger	Translate Spreadsheet Cell Ranges to Rows and Columns
CellularAutomaton	One-Dimensional Cellular Automata
cellVolumeDist	Functions to fit cell volume distributions and thereby estimate cell growth rates and division times
cem	Coarsened Exact Matching
cems	Conditional Expectation Manifolds
CensMixReg	Censored Linear Mixture Regression Models
censNID	censored NID samples
censorcopula	Estimate Parameter of Bivariate Copula
censReg	Censored Regression (Tobit) Models
CensRegMod	Fits Normal and Student-t Censored Regression Model
censusr	Collect Data from the Census API
cents	Censored time series
CEoptim	Cross-Entropy R Package for Optimization
CePa	Centrality-based pathway enrichment
CepLDA	Discriminant Analysis of Time Series in the Presence of Within-Group Spectral Variability
cepp	Context Driven Exploratory Projection Pursuit
CerioliOutlierDetection	Outlier detection using the iterated RMCD method of Cerioli (2010)
cernn	Covariance Estimation Regularized by Nuclear Norm Penalties
CFC	Cause-Specific Framework for Competing-Risk Analysis
CfEstimateQuantiles	Estimate quantiles using any order Cornish-Fisher expansion
cffdrs	Canadian Forest Fire Danger Rating System
cg	Compare Groups, Analytically and Graphically
cgam	Constrained Generalized Additive Model
cgAUC	Calculate AUC-type measure when gold standard is continuous and the corresponding optimal linear combination of variables with respect to it
cgdsr	R-Based API for Accessing the MSKCC Cancer Genomics Data Server (CGDS)
cggd	Continuous Generalized Gradient Descent
cgh	Microarray CGH analysis using the Smith-Waterman algorithm
cghFLasso	Detecting hot spot on CGH array data with fused lasso regression
cghseg	Segmentation Methods for Array CGH Analysis
CGP	Composite Gaussian process models
cgwtools	Miscellaneous Tools
ChainLadder	Statistical Methods and Models for Claims Reserving in General Insurance
changepoint	Methods for Changepoint Detection
ChannelAttribution	Markov Model for the Online Multi-Channel Attribution Problem
ChannelAttributionApp	Shiny Web Application for the Multichannel Attribution Problem
ChaosGame	Chaos Game
ChargeTransport	Charge Transfer Rates and Charge Carrier Mobilities
CHAT	Clonal Heterogeneity Analysis Tool
CHCN	Canadian Historical Climate Network
cheb	Discrete Linear Chebyshev Approximation
chebpol	Multivariate Chebyshev Interpolation
CheckDigit	Calculate and verify check digits
checkmate	Fast and Versatile Argument Checks
checkpoint	Install Packages from Snapshots on the Checkpoint Server for Reproducibility
cheddar	Analysis and Visualisation of Ecological Communities
chemCal	Calibration Functions for Analytical Chemistry
chemometrics	Multivariate Statistical Analysis in Chemometrics
ChemometricsWithR	Chemometrics with R - Multivariate Data Analysis in the Natural Sciences and Life Sciences
ChemometricsWithRData	Data for package ChemometricsWithR
ChemoSpec	Exploratory Chemometrics for Spectroscopy
cherry	Multiple Testing Methods for Exploratory Research
childsds	Calculation of standard deviation scores adduced from different growth standards
chillR	Statistical Methods for Phenology Analysis in Temperate Fruit Trees
chipPCR	Toolkit of Helper Functions to Pre-Process Amplification Data
chngpt	Change Point Regression
CHNOSZ	Chemical Thermodynamics and Activity Diagrams
choiceDes	Design Functions for Choice Studies
ChoiceModelR	Choice Modeling in R
choplump	Choplump tests
chopthin	The Chopthin Resampler
chords	Estimation in respondent driven samples
choroplethr	Simplify the Creation of Choropleth Maps in R
choroplethrAdmin1	Contains an Administrative-Level-1 Map of the World
choroplethrMaps	Contains maps used by the choroplethr package
chromer	Interface to Chromosome Counts Database API
chromoR	Analysis of chromosomal interactions data (correction, segmentation and comparison)
chron	Chronological Objects which can Handle Dates and Times
CHsharp	Choi and Hall Style Data Sharpening
chunked	Chunkwise Text-File Processing for 'dplyr'
CIDnetworks	Generative Models for Complex Networks with Conditionally Independent Dyadic Structure
CIFsmry	Weighted summary of cumulative incidence functions
cin	Causal Inference for Neuroscience
CINID	Curculionidae INstar IDentification
CINOEDV	Co-Information based N-Order Epistasis Detector and Visualizer
CircE	Circumplex models Estimation
circlize	Circular Visualization
CircNNTSR	CircNNTSR: An R package for the statistical analysis of circular data using nonnegative trigonometric sums (NNTS) models
CircOutlier	Detection of Outliers in Circular-Circular Regression
CircStats	Circular Statistics, from "Topics in circular Statistics" (2001)
circular	Circular Statistics
cIRT	Choice Item Response Theory
cit	Causal Inference Test
CITAN	CITation ANalysis Toolpack
citbcmst	CIT Breast Cancer Molecular SubTypes Prediction
citccmst	CIT Colon Cancer Molecular SubTypes Prediction
CityPlot	Visualization of structure and contents of a database
cjoint	AMCE Estimator for Conjoint Experiments
ckanr	Client for the Comprehensive Knowledge Archive Network ('CKAN') 'API'
Ckmeans.1d.dp	Optimal k-Means Clustering for One-Dimensional Data
cladoRcpp	C++ implementations of phylogenetic cladogenesis calculations
ClamR	Time Series Modeling for Climate Change Proxies
clarifai	Access to Clarifai API
class	Functions for Classification
classGraph	Construct Graphs of S4 Class Hierarchies
classifly	Explore classification models in high dimensions
classify	Classification Accuracy and Consistency under IRT models
classInt	Choose Univariate Class Intervals
classyfire	Robust multivariate classification using highly optimised SVM ensembles
cleangeo	Cleaning Geometries from Spatial Objects
clere	Simultaneous Variables Clustering and Regression
clhs	Conditioned Latin Hypercube Sampling
ClickClust	Model-Based Clustering of Categorical Sequences
clickstream	Analyzes Clickstreams Based on Markov Chains
clifro	Easily Download and Visualise Climate Data from CliFlo
ClimClass	Climate Classification According to Several Indices
climdex.pcic	PCIC Implementation of Climdex Routines
clime	Constrained L1-minimization for Inverse (covariance) Matrix Estimation
climtrends	Statistical Methods for Climate Sciences
climwin	Climate Window Analysis
clinfun	Clinical Trial Design and Data Analysis Functions
clinsig	Clinical Significance Functions
clinUtiDNA	Clinical Utility of DNA Testing
clipr	Read and Write from the System Clipboard
clisymbols	Unicode Symbols at the R Prompt
CLME	Constrained Inference for Linear Mixed Effects Models
clogitboost	Boosting Conditional Logit Model
clogitL1	Fitting exact conditional logistic regression with lasso and elastic net penalties
cloudUtil	Cloud Util Plots
clpAPI	R Interface to C API of COIN-OR Clp
CLSOCP	A smoothing Newton method SOCP solver
clttools	Central Limit Theorem Experiments (Theoretical and Simulation)
clue	Cluster Ensembles
ClueR	CLUster Evaluation (CLUE)
clues	Clustering Method Based on Local Shrinking
CluMix	Clustering and Visualization of Mixed-Type Data
clusrank	Wilcoxon Rank Sum Test for Clustered Data
cluster	"Finding Groups in Data": Cluster Analysis Extended Rousseeuw et al.
clusterCrit	Clustering Indices
cluster.datasets	Cluster Analysis Data Sets
clusterfly	Explore clustering interactively using R and GGobi
clusterGeneration	Random Cluster Generation (with Specified Degree of Separation)
clusterGenomics	Identifying clusters in genomics data by recursive partitioning
clustering.sc.dp	Optimal Distance-Based Clustering for Multidimensional Data with Sequential Constraint
clusterPower	Power calculations for cluster-randomized and cluster-randomized crossover trials
clusterRepro	Reproducibility of gene expression clusters
clusterSEs	Calculate Cluster-Robust p-Values and Confidence Intervals
clusterSim	Searching for Optimal Clustering Procedure for a Data Set
ClusterStability	Assessment of Stability of Individual Objects or Clusters in Partitioning Solutions
clustertend	Check the Clustering Tendency
clusteval	Evaluation of Clustering Algorithms
ClustGeo	Clustering of Observations with Geographical Constraints
clustMD	Model Based Clustering for Mixed Data
clustMixType	k-Prototypes Clustering for Mixed Variable-Type Data
ClustMMDD	Variable Selection in Clustering by Mixture Models for Discrete Data
ClustOfVar	Clustering of variables
clustrd	Methods for joint dimension reduction and clustering
clustsig	Significant Cluster Analysis
ClustVarLV	Clustering of Variables Around Latent Variables
clustvarsel	Variable Selection for Model-Based Clustering
clv	Cluster Validation Techniques
clValid	Validation of Clustering Results
cmaes	Covariance Matrix Adapting Evolutionary Strategy
cmaesr	Covariance Matrix Adaptation Evolution Strategy
CMC	Cronbach-Mesbah Curve
CMF	Collective matrix factorization
cmm	Categorical Marginal Models
cmna	Computational Methods for Numerical Analysis
CMPControl	Control Charts for Conway-Maxwell-Poisson Distribution
CMplot	Circle Manhattan Plot
cmprsk	Subdistribution Analysis of Competing Risks
cmprskQR	Analysis of Competing Risks Using Quantile Regressions
cmrutils	Misc Functions of the Center for the Mathematical Research
cmsaf	Tools for CM SAF Netcdf Data
cmvnorm	The Complex Multivariate Gaussian Distribution
cna	A Package for Coincidence Analysis (CNA)
cncaGUI	Canonical Non-Symmetrical Correspondence Analysis in R
cnmlcd	Maximum Likelihood Estimation of a Log-Concave Density Function
CNOGpro	Copy Numbers of Genes in prokaryotes
CNprep	Pre-process DNA Copy Number (CN) Data for Detection of CN Events
CNVassoc	Association Analysis of CNV Data and Imputed SNPs
CNVassocData	Example data sets for association analysis of CNV data
coala	A Framework for Coalescent Simulation
coalescentMCMC	MCMC Algorithms for the Coalescent
coarseDataTools	A Collection of Functions to Help with Analysis of Coarsely Observed Data
COBRA	Nonlinear Aggregation of Predictors
cobs	COnstrained B-Splines (Sparse Matrix Based)
CoClust	Copula Based Cluster Analysis
cocor	Comparing Correlations
cocorresp	Co-Correspondence Analysis Methods
cocron	Statistical Comparisons of Two or more Alpha Coefficients
coda	Output Analysis and Diagnostics for MCMC
codadiags	Markov chain Monte Carlo burn-in based on "bridge" statistics
cOde	Automated C Code Generation for Use with the 'deSolve' and 'bvpSolve' Packages
codep	Multiscale Codependence Analysis
codetools	Code Analysis Tools for R
codingMatrices	Alternative Factor Coding Matrices for Linear Model Formulae
codyn	Community Dynamics Metrics
coefficientalpha	Robust Coefficient Alpha and Omega with Missing and Non-Normal Data
coefplot	Plots Coefficients from Fitted Models
coenocliner	Coenocline Simulation
coenoflex	Gradient-Based Coenospace Vegetation Simulator
coexist	Species coexistence modeling and analysis
cofeatureR	Generate Cofeature Matrices
CoImp	Copula based imputation method
coin	Conditional Inference Procedures in a Permutation Test Framework
CoinMinD	Simultaneous Confidence Interval for Multinomial Proportion
cold	Count Longitudinal Data
CollocInfer	Collocation Inference for Dynamic Systems
coloc	Colocalisation tests of two genetic traits
coloredICA	Implementation of Colored Independent Component Analysis and Spatial Colored Independent Component Analysis
colorfulVennPlot	Plot and add custom coloring to Venn diagrams for 2-dimensional, 3-dimensional and 4-dimensional data
colorhcplot	Colorful Hierarchical Clustering Dendrograms
ColorPalette	Color Palettes Generator
colorRamps	Builds color tables
colorscience	Color Science Methods and Data
colorspace	Color Space Manipulation
colortools	Tools for colors in a Hue-Saturation-Value (HSV) color model
colourlovers	R Client for the COLOURlovers API
comato	Analysis of Concept Maps
COMBIA	Synergy/Antagonism Analyses of Drug Combinations
combinat	combinatorics utilities
Combine	Game-Theoretic Probability Combination
CombinePValue	Combine a Vector of Correlated p-values
CombinS	Construction Methods of some Series of PBIB Designs
CombMSC	Combined Model Selection Criteria
comclim	Community climate statistics
cometExactTest	Exact Test from the Combinations of Mutually Exclusive Alterations (CoMEt) Algorithm
comf	Calculation of Thermal Comfort Indices
ComICS	Computational Methods for Immune Cell-Type Subsets
commandr	Command pattern in R
commentr	Print Nicely Formatted Comments for Use in Script Files
CommonJavaJars	Useful libraries for building a Java based GUI under R
commonmark	Bindings to the CommonMark Reference Implementation
CommonTrend	Extract and plot common trends from a cointegration system. Calculate P-value for Johansen Statistics
CommT	Comparative Phylogeographic Analysis using the Community Tree Framework
COMMUNAL	Robust Selection of Cluster Number K
CommunityCorrelogram	Ecological Community Correlogram
Comp2ROC	Compare Two ROC Curves that Intersect
compactr	Creates empty plots with compact axis notation
compare	Comparing Objects for Differences
compareC	Compare Two Correlated C Indices with Right-censored Survival Outcome
CompareCausalNetworks	Interface to Diverse Estimation Methods of Causal Networks
compareDF	Do a Git Style Diff of the Rows Between Two Dataframes with Similar Structure
compareGroups	Descriptive Analysis by Groups
compareODM	comparison of medical forms in CDISC ODM format
CompareTests	Correct for Verification Bias in Diagnostic Accuracy & Agreement
comparison	Multivariate likelihood ratio calculation and evaluation
compeir	Event-specific incidence rates for competing risks data
compendiumdb	Tools for Retrieval and Storage of Functional Genomics Data
CompGLM	Conway-Maxwell-Poisson GLM and distribution functions
compHclust	Complementary Hierarchical Clustering
Compind	Composite Indicators Functions
ComplexAnalysis	Numerically Evaluate Integrals and Derivatives (also Higher Order) of Vector- And Complex-Valued Functions
complmrob	Robust Linear Regression with Compositional Data as Covariates
CompLognormal	Functions for actuarial scientists
compoisson	Conway-Maxwell-Poisson Distribution
COMPoissonReg	Conway-Maxwell Poisson (COM-Poisson) Regression
Compositional	Compositional Data Analysis
compositions	Compositional Data Analysis
compound.Cox	Regression Estimation Based on the Compound Covariate Method Under the Cox Proportional Hazard Model
Compounding	Computing Continuous Distributions
CompQuadForm	Distribution function of quadratic forms in normal variables
CompR	Paired Comparison Data Analysis
CompRandFld	Composite-Likelihood Based Analysis of Random Fields
compute.es	Compute Effect Sizes
Conake	Continuous Associated Kernel Estimation
ConConPiWiFun	Optimisation with Continuous Convex Piecewise (Linear and Quadratic) Functions
concor	Concordance
concordance	Product Concordance
concreg	Concordance regression
cond	Approximate conditional inference for logistic and loglinear models
condformat	Conditional Formatting in Data Frames
condGEE	Parameter estimation in conditional GEE for recurrent event gap times
condmixt	Conditional Density Estimation with Neural Network Conditional Mixtures
condMVNorm	Conditional Multivariate Normal Distribution
CONDOP	Condition-Dependent Operon Predictions
CondReg	Condition Number Regularized Covariance Estimation
condvis	Conditional Visualization for Statistical Models
coneproj	Primal or Dual Cone Projections with Routines for Constrained Regression
conf.design	Construction of factorial designs
confidence	Confidence Estimation of Environmental State Classifications
conformal	Conformal Prediction for Regression and Classification
confreq	Configural Frequencies Analysis Using Log-Linear Modeling
conicfit	Algorithms for Fitting Circles, Ellipses and Conics Based on the Work by Prof. Nikolai Chernov
conics	Plot Conics
conjoint	Conjoint analysis package
ConjointChecks	A package to check the cancellation axioms of conjoint measurement
connect3	A Tool for Reproducible Research by Converting 'LaTeX' Files Generated by R Sweave to Rich Text Format Files
ConnMatTools	Tools for working with connectivity matrices
conover.test	Conover-Iman Test of Multiple Comparisons Using Rank Sums
ConSpline	Partial Linear Least-Squares Regression using Constrained Splines
ConsRank	Compute the Median Ranking(s) According to the Kemeny's Axiomatic Approach
constrainedKriging	Constrained, Covariance-Matching Constrained and Universal Point or Block Kriging
ContaminatedMixt	Model-Based Clustering and Classification with the Multivariate Contaminated Normal Distribution
contfrac	Continued fractions
conting	Bayesian Analysis of Contingency Tables
contoureR	Contouring of Non-Regular Three-Dimensional Data
contrast	A collection of contrast methods
controlTest	Median Comparison for Two-Sample Right-Censored Survival Data
ConvCalendar	Converts dates between calendars
ConvergenceConcepts	Seeing convergence concepts in action
convevol	Quantifies and assesses the significance of convergent evolution
convoSPAT	Convolution-Based Nonstationary Spatial Modeling
cooccur	Probabilistic Species Co-Occurrence Analysis in R
coop	Co-Operation: Fast Covariance, Correlation, and Cosine Similarity Operations
cooptrees	Cooperative aspects of optimal trees in weighted graphs
copBasic	General Bivariate Copula Theory and Many Utility Functions
copCAR	Fitting the copCAR Regression Model for Discrete Areal Data
cope	Coverage Probability Excursion (CoPE) Sets
copula	Multivariate Dependence with Copulas
CopulaDTA	Copula Based Bivariate Beta-Binomial Model for Diagnostic Test Accuracy Studies
copulaedas	Estimation of Distribution Algorithms Based on Copulas
Copula.Markov	Estimation and Statistical Process Control Under Copula-Based Time Series Models
CopulaRegression	Bivariate Copula Based Regression Models
CopulaREMADA	Copula Mixed Effect Models for Bivariate and Trivariate Meta-Analysis of Diagnostic Test Accuracy Studies
CopyDetect	Computing Statistical Indices to Detect Answer Copying on Multiple-Choice Tests
CopyNumber450kCancer	Baseline Correction for Copy Number Data from Cancer Samples
corclass	Correlational Class Analysis
corcounts	Generate correlated count random variables
cord	Community Estimation in G-Models via CORD
CORE	Cores of Recurrent Events
CORElearn	Classification, Regression and Feature Evaluation
coreNLP	Wrappers Around Stanford CoreNLP Tools
coreTDT	TDT for compound heterozygous and recessive models
corHMM	Analysis of Binary Character Evolution
corkscrew	Preprocessor for Data Modeling
CORM	The Clustering of Regression Models Method
corpcor	Efficient Estimation of Covariance and (Partial) Correlation
corpora	Statistics and data sets for corpus frequency data
CorrBin	Nonparametrics with Clustered Binary and Multinomial Data
CorReg	Linear Regression Based on Linear Structure Between Covariates
corregp	Functions and Methods for Correspondence Regression
Correlplot	A collection of functions for graphing correlation matrices
corrgram	Plot a Correlogram
CorrMixed	Estimate Correlations Between Repeatedly Measured Endpoints (E.g., Reliability) Based on Linear Mixed-Effects Models
corrplot	Visualization of a correlation matrix
corrsieve	CorrSieve
corTools	Tools for processing data after a Genome Wide Association Study
COSINE	COndition SpecIfic sub-NEtwork
cosinor	Tools for estimating and predicting the cosinor model
cosmoFns	Functions for cosmological distances, times, luminosities, etc
CosmoPhotoz	Photometric redshift estimation using generalized linear models
cosmosR	cosmosR
cosso	Fit Regularized Nonparametric Regression Models Using COSSO Penalty
costat	Time series costationarity determination
cotrend	Consistant Cotrend Rank Selection
couchDB	Connect and work with couchDB databases
COUNT	Functions, data and code for count data
Countr	Flexible Univariate and Bivariate Count Process Probability
countrycode	Convert Country Names and Country Codes
CountsEPPM	Mean and Variance Modeling of Count Data
COUSCOus	A Residue-Residue Contact Detecting Method
covafillr	Local Polynomial Regression of State Dependent Covariates in State-Space Models
covBM	Brownian Motion Processes for 'nlme'-Models
covLCA	Latent Class Models with Covariate Effects on Underlying and Measured Variables
covmat	Covariance Matrix Estimation
covr	Test Coverage for Packages
covreg	A simultaneous regression model for the mean and covariance
covRobust	Robust Covariance Estimation via Nearest Neighbor Cleaning
CovSel	Model-Free Covariate Selection
covsep	Tests for Determining if the Covariance Structure of 2-Dimensional Data is Separable
covTest	Computes covariance test for adaptive linear modelling
cowplot	Streamlined Plot Theme and Plot Annotations for 'ggplot2'
cowsay	Messages, Warnings, Strings with Ascii Animals
CoxBoost	Cox models by likelihood based boosting for a single survival endpoint or competing risks
coxinterval	Cox-Type Models for Interval-Censored Data
coxme	Mixed Effects Cox Models
Coxnet	Regularized Cox Model
coxphf	Cox regression with Firth's penalized likelihood
coxphw	Weighted Estimation in Cox Regression
CoxPlus	Cox Regression (Proportional Hazards Model) with Multiple Causes and Mixed Effects
CoxRidge	Cox Models with Dynamic Ridge Penalties
coxrobust	Robust Estimation in Cox Model
coxsei	Fitting a CoxSEI Model
CP	Conditional Power Calculations
cp4p	Calibration Plot for Proteomics
cpa	Confirmatory Path Analysis through the d-sep tests
cpca	Methods to perform Common Principal Component Analysis (CPCA)
CPE	Concordance Probability Estimates in Survival Analysis
CpGassoc	Association Between Methylation and a Phenotype of Interest
cpgen	Parallelized Genomic Prediction and GWAS
CpGFilter	CpG Filtering Method Based on Intra-class Correlation Coefficients
CPHshape	Find the maximum likelihood estimator of the shape constrained hazard baseline and the effect parameters in the Cox proportional hazards model
cpk	Clinical Pharmacokinetics
cplexAPI	R Interface to C API of IBM ILOG CPLEX
cplm	Compound Poisson Linear Models
cpm	Sequential and Batch Change Detection Using Parametric and Nonparametric Methods
CPMCGLM	Correction of the pvalue after multiple coding
Cprob	Conditional probability function of a competing event
cqrReg	Quantile, Composite Quantile Regression and Regularized Versions
cquad	Conditional Maximum Likelihood for Quadratic Exponential Models for Binary Panel Data
CR	Power Calculation for Weighted Log-Rank Tests in Cure Rate Models
CRAC	Cosmology R Analysis Code
crackR	Probabilistic damage tolerance analysis for fatigue cracking of metallic aerospace structures
cramer	Multivariate nonparametric Cramer-Test for the two-sample-problem
crandatapkgs	Find Data-Only Packages on CRAN
crank	Completing Ranks
cranlogs	Download Logs from the 'RStudio' 'CRAN' Mirror
crantastic	Various R tools for http://crantastic.org/
crawl	Fit Continuous-Time Correlated Random Walk Models to Animal Movement Data
crayon	Colored Terminal Output
crblocks	Categorical Randomized Block Data Analysis
crch	Censored Regression with Conditional Heteroscedasticity
CreditMetrics	Functions for calculating the CreditMetrics risk model
creditr	Credit Default Swaps in R
credule	Credit Default Swap Functions
CRF	CRF - Conditional Random Fields
cricketr	Analyze Cricketers Based on ESPN Cricinfo Statsguru
crimCV	Group-Based Modelling of Longitudinal Data
crimelinkage	Statistical Methods for Crime Series Linkage
CRM	Continual Reassessment Method (CRM) for Phase I Clinical Trials
crmn	CCMN and other noRMalizatioN methods for metabolomics data
crmPack	Object-Oriented Implementation of CRM Designs
crn	Downloads and Builds datasets for Climate Reference Network
crop	Graphics Cropping Tool
crossdes	Construction of Crossover Designs
crossmatch	The Cross-match Test
Crossover	Analysis and Search of Crossover Designs
crossReg	Confidence intervals for crossover points of two simple regression lines
crossval	Generic Functions for Cross Validation
crp.CSFP	CreditRisk+ Portfolio Model
crqa	Cross-Recurrence Quantification Analysis for Categorical and Continuous Time-Series
crrp	Penalized Variable Selection in Competing Risks Regression
crrSC	Competing risks regression for Stratified and Clustered data
crrstep	Stepwise Covariate Selection for the Fine & Gray Competing Risks Regression Model
crs	Categorical Regression Splines
crskdiag	Diagnostics for Fine and Gray Model
CRTgeeDR	Doubly Robust Inverse Probability Weighted Augmented GEE Estimator
CRTSize	Sample Size Estimation Functions for Cluster Randomized Trials
crunch	Crunch.io Data Tools
cruts	Interface to Climatic Research Unit Time-Series Version 3.21 Data
CrypticIBDcheck	Identifying cryptic relatedness in genetic association studies
CryptRndTest	Statistical Tests for Cryptographic Randomness
csampling	Functions for Conditional Simulation in Regression-Scale Models
cSFM	Covariate-adjusted Skewed Functional Model (cSFM)
cshapes	CShapes Dataset and Utilities
cslogistic	Conditionally Specified Logistic Regression
csn	Closed Skew-Normal Distribution
csrplus	Methods to Test Hypotheses on the Distribution of Spatial Point Processes
csSAM	csSAM - cell-specific Significance Analysis of Microarrays
cstar	Substantive significance testing for regression estimates and marginal effects
csvread	Fast Specialized CSV File Loader
ctmcmove	Modeling Animal Movement with Continuous-Time Discrete-Space Markov Chains
ctmm	Continuous-Time Movement Modeling
cts	Continuous Time Autoregressive Models
ctsem	Continuous Time Structural Equation Modelling
CTT	Classical Test Theory Functions
CTTShiny	Classical Test Theory via Shiny
ctv	CRAN Task Views
CUB	A Class of Mixture Models for Ordinal Data
cubature	Adaptive multivariate integration over hypercubes
cubfits	Codon Usage Bias Fits
Cubist	Rule- and Instance-Based Regression Modeling
cudaBayesreg	CUDA Parallel Implementation of a Bayesian Multilevel Model for fMRI Data Analysis
cudaBayesregData	Data sets for the examples used in the package "cudaBayesreg"
cudia	CUDIA Cross-level Imputation
CUFF	Charles's Utility Function using Formula
CUMP	Analyze Multivariate Phenotypes by Combining Univariate results
cumplyr	Extends ddply to allow calculation of cumulative quantities
cumSeg	Change point detection in genomic sequences
curl	A Modern and Flexible Web Client for R
currentSurvival	Estimation of CCI and CLFS Functions
curvetest	The package will formally test two curves represented by discrete data sets to be statistically equal or not when the errors of the two curves were assumed either equal or not using the tube formula to calculate the tail probabilities
curvHDR	Filtering of Flow Cytometry Samples
cusp	Cusp-Catastrophe Model Fitting Using Maximum Likelihood
customizedTraining	Customized Training for Lasso and Elastic-Net Regularized Generalized Linear Models
CUSUMdesign	Compute Decision Interval and Average Run Length for CUSUM Charts
cutoffR	CUTOFF: A Spatio-temporal Imputation Method
cuttlefish.model	An R package to perform LPUE standardization and stock assessment of the English Channel cuttlefish stock using a two-stage biomass model
cvAUC	Cross-Validated Area Under the ROC Curve Confidence Intervals
CVcalibration	Estimation of the Calibration Equation with Error-in Observations
cvplogistic	Penalized Logistic Regression Model using Majorization Minimization by Coordinate Descent (MMCD) Algorithm
cvq2	Calculate the predictive squared correlation coefficient
CVST	Fast Cross-Validation via Sequential Testing
CVThresh	Level-Dependent Cross-Validation Thresholding
cvTools	Cross-validation tools for regression models
CVTuningCov	Regularized Estimators of Covariance Matrices with CV Tuning
cvxbiclustr	Convex Biclustering Algorithm
cvxclustr	Splitting methods for convex clustering
cwhmisc	Miscellaneous Functions for Math, Plotting, Printing, Statistics, Strings, and Tools
cwm	Cluster Weighted Models by EM algorithm
cxxfunplus	extend cxxfunction by saving the dynamic shared objects
cycleRtools	Tools for Cycling Data Analysis
cyclocomp	Cyclomatic Complexity of R Code
cycloids	cycloids
cymruservices	Query 'Team Cymru' 'IP' Address, Autonomous System Number ('ASN'), Border Gateway Protocol ('BGP'), Bogon and 'Malware' Hash Data Services
cyphid	Cycle and Phase Identification for mastication data
cytoDiv	Cytometric diversity indices
D2C	Predicting Causal Direction from Dependency Features
d3heatmap	Interactive Heat Maps Using 'htmlwidgets' and 'D3.js'
D3M	Two Sample Test with Wasserstein Metric
d3Network	Tools for creating D3 JavaScript network, tree, dendrogram, and Sankey graphs from R
DAAG	Data Analysis and Graphics Data and Functions
DAAGbio	Data Sets and Functions, for demonstrations with expression arrays and gene sequences
DAAGxtras	Data Sets and Functions, supplementary to DAAG
dad	Three-Way Data Analysis Through Densities
dae	Functions Useful in the Design and ANOVA of Experiments
daewr	Design and Analysis of Experiments with R
daff	Diff, Patch and Merge for Data.frames
dafs	Data analysis for forensic scientists
dagbag	Learning directed acyclic graphs (DAGs) through bootstrap aggregating
DAGGER	Consensus genetic maps
dagitty	Graphical Analysis of Structural Causal Models
dagR	R functions for directed acyclic graphs
Daim	Diagnostic accuracy of classification models
DAISIE	Dynamical Assembly of Islands by Speciation, Immigration and Extinction
DAKS	Data Analysis and Knowledge Spaces
DALY	DALY Calculator - A GUI for stochastic DALY calculation in R
DAMisc	Dave Armstrong's Miscellaneous Functions
DAMOCLES	Dynamic Assembly Model of Colonization, Local Extinction and Speciation
dams	Dams in the United States from the National Inventory of Dams (NID)
DandEFA	Dandelion Plot for R-Mode Exploratory Factor Analysis
darch	Package for Deep Architectures and Restricted Boltzmann Machines
Dark	The Analysis of Dark Adaptation Data
darts	Statistical Tools to Analyze Your Darts Game
dashboard	Interactive Data Visualization with D3.js
DatABEL	File-Based Access to Large Matrices Stored on HDD in Binary Format
datacheck	Tools for Checking Data Consistency
datacheckr	Data Checking
DataClean	Data Cleaning
DataCombine	Tools for Easily Combining and Cleaning Data Sets
datadr	Divide and Recombine for Large, Complex Data
DataExplorer	Data Explorer
dataframes2xls	dataframes2xls writes data frames to xls files
datafsm	Estimating Finite State Machine Models from Data
DataLoader	Import Multiple File Types
datamap	A system for mapping foreign objects to R variables and environments
datamart	Unified access to your data sources
dataonderivatives	Easily Source Publicly Available Data on Derivatives
datapack	A Flexible Container to Transport and Manipulate Data and Associated Resources
dataQualityR	Performs variable level data quality checks and generates summary statistics
dataRetrieval	Retrieval Functions for USGS and EPA Hydrologic and Water Quality Data
datastepr	An Implementation of a SAS-Style Data Step
data.table	Extension of Data.frame
data.tree	General Purpose Hierarchical Data Structure
datautils	Support functions for packages VBmix, semisupKernelPCA, and patchPlot
dataview	Data and Workspace Browser for Terminals
date	Functions for handling dates
DATforDCEMRI	Deconvolution Analysis Tool for Dynamic Contrast Enhanced MRI
dave	Functions for "Data Analysis in Vegetation Ecology"
Davies	The Davies quantile function
dawai	Discriminant Analysis with Additional Information
dbarts	Discrete Bayesian Additive Regression Trees Sampler
dbConnect	Provides a graphical user interface to connect with databases that use MySQL
dbEmpLikeGOF	Goodness-of-fit and two sample comparison tests using sample entropy
dbEmpLikeNorm	Test for joint assessment of normality
DBGSA	methods of distance-based gene set functional enrichment analysis
DBI	R Database Interface
DBKGrad	Discrete Beta Kernel Graduation of Mortality Data
dblcens	Compute the NPMLE of distribution from doubly censored data
dbmss	Distance-Based Measures of Spatial Structures
dbscan	Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms
dbstats	Distance-Based Statistics
DCchoice	Analyzing Dichotomous Choice Contingent Valuation Data
dcemriS4	A Package for Image Analysis of DCE-MRI (S4 Implementation)
DCGL	Differential Co-expression Analysis and Differential Regulation Analysis of Gene Expression Microarray Data
dcGOR	Analysis of Ontologies and Protein Domain Annotations
dChipIO	Methods for Reading dChip Files
DCL	Claims Reserving under the Double Chain Ladder Model
dclone	Data Cloning and MCMC Tools for Maximum Likelihood Methods
DCluster	Functions for the Detection of Spatial Clusters of Diseases
dcmle	Hierarchical Models Made Easy with Data Cloning
dcmr	Attribute profile estimation using Diagnostic Classification Models and MCMC
DCODE	List Linear n-Peptide Constraints for Overlapping Protein Regions
dCovTS	Distance Covariance and Correlation for Time Series Analysis
dcv	Conventional Cross-validation statistics for climate-growth model
ddalpha	Depth-Based Classification and Calculation of Data Depth
DDD	Diversity-Dependent Diversification
ddeploy	Wrapper for the Duke Deploy REST API
DDHFm	Variance Stabilization by Data-Driven Haar-Fisz (for Microarrays)
DDIwR	DDI with R
ddpcr	Analysis and Visualization of Droplet Digital PCR in R and on the Web
ddR	Distributed Data Structures in R
DDRTree	Learning Principal Graphs with DDRTree
ddst	Data driven smooth test
deal	Learning Bayesian Networks with Mixed Variables
deamer	Deconvolution density estimation with adaptive methods for a variable prone to measurement error
debug	MVB's debugger for R
DECIDE	DEComposition of Indirect and Direct Effects
DecisionCurve	Calculate and Plot Decision Curves
decisionSupport	Quantitative Support of Decision Making under Uncertainty
decode	Differential Co-Expression and Differential Expression Analysis
decompr	Global-Value-Chain Decomposition
decon	Deconvolution Estimation in Measurement Error Models
deconstructSigs	Identifies Signatures Present in a Tumor Sample
Deducer	A Data Analysis GUI for R
DeducerExtras	Additional dialogs and functions for Deducer
DeducerPlugInExample	Deducer Plug-in Example
DeducerPlugInScaling	Reliability and factor analysis plugin
DeducerSpatial	Deducer for spatial data analysis
DeducerSurvival	Add Survival Dialogue to Deducer
DeducerText	Deducer GUI for Text Data
deducorrect	Deductive Correction, Deductive Imputation, and Deterministic Correction
deductive	Data Correction and Imputation Using Deductive Methods
deepboost	Deep Boosting Ensemble Modeling
deepnet	deep learning toolkit in R
DEEPR	Dirichlet-multinomial Evolutionary Event Profile Randomization (DEEPR) test
deformula	Integration of One-Dimensional Functions with Double Exponential Formulas
degenes	Detection of differentially expressed genes
degreenet	Models for Skewed Count Distributions Relevant to Networks
Delaporte	Statistical Functions for the Delaporte Distribution
deldir	Delaunay Triangulation and Dirichlet (Voronoi) Tessellation
delt	Estimation of Multivariate Densities Using Adaptive Partitions
deltaPlotR	Identification of dichotomous differential item functioning (DIF) using Angoff's Delta Plot method
Demerelate	Functions to calculate relatedness on diploid genetic data
DEMEtics	Evaluating the genetic differentiation between populations based on Gst and D values
demi	Differential Expression from Multiple Indicators
deming	Deming, Thiel-Sen and Passing-Bablock Regression
demography	Forecasting mortality, fertility, migration and population data
demoKde	Kernel Density Estimation for Demonstration Purposes
DEMOVA	DEvelopment (of Multi-Linear QSPR/QSAR) MOdels VAlidated using Test Set
dendextend	Extending R's Dendrogram Functionality
dendextendRcpp	Faster Dendrogram Manipulation using 'Rcpp'
dendroextras	Extra functions to cut, label and colour dendrogram clusters
dendrometeR	Analyzing Dendrometer Data
DendSer	Dendrogram seriation: ordering for visualisation
dendsort	Modular Leaf Ordering Methods for Dendrogram Nodes
denovolyzeR	Statistical Analyses of De Novo Genetic Variants
denpro	Visualization of Multivariate Functions, Sets, and Data
densityClust	Clustering by Fast Search and Find of Density Peaks
Density.T.HoldOut	Density.T.HoldOut: Non-combinatorial T-estimation Hold-Out for density estimation
denstrip	Density strips and other methods for compactly illustrating distributions
DEoptim	Global Optimization by Differential Evolution
DEoptimR	Differential Evolution Optimization in Pure R
depend.truncation	Statistical Inference for Parametric and Semiparametric Models Based on Dependently Truncated Data
depmix	Dependent Mixture Models
depmixS4	Dependent Mixture Models - Hidden Markov Models of GLMs and Other Distributions in S4
depth	Depth functions tools for multivariate analysis
depth.plot	Multivariate Analogy of Quantiles
DepthProc	Statistical Depth Functions for Multivariate Analysis
depthTools	Depth Tools Package
dequer	An R 'Deque' Container
Deriv	Symbolic Differentiation
descomponer	Seasonal Adjustment by Frequency Analysis
descr	Descriptive Statistics
DescribeDisplay	An Interface to the DescribeDisplay GGobi Plugin
describer	Describe Data in R Using Common Descriptive Statistics
DescTools	Tools for Descriptive Statistics
deseasonalize	Optimal deseasonalization for geophysical time series using AR fitting
designGG	Computational tool for designing genetical genomics experiments
designGLMM	Finding Optimal Block Designs for a Generalised Linear Mixed Model
designmatch	Construction of Optimally Matched Samples for Randomized Experiments and Observational Studies that are Balanced by Design
desiR	Desirability Functions for Ranking, Selecting, and Integrating Data
desirability	Desirability Function Optimization and Ranking
desire	Desirability functions in R
DESnowball	Bagging with Distance-based Regression for Differential Gene Expression Analyses
deSolve	Solvers for Initial Value Problems of Differential Equations (ODE, DAE, DDE)
DESP	Estimation of Diagonal Elements of Sparse Precision-Matrices
desplot	Plotting Field Plans for Agricultural Experiments
detect	Analyzing Wildlife Data with Detection Error
detector	Detect Data Containing Personally Identifiable Information
deTestSet	Testset for differential equations
DetMCD	Implementation of the DetMCD Algorithm (Robust and Deterministic Estimation of Location and Scatter)
DetR	Suite of Deterministic and Robust Algorithms for Linear Regression
detrendeR	Start the detrendeR Graphical User Interface (GUI)
DetSel	A computer program to detect markers responding to selection
devEMF	EMF Graphics Output Device
Devore7	Data sets from Devore's "Prob and Stat for Eng (7th ed)"
devtools	Tools to Make Developing R Packages Easier
df2json	Convert a dataframe to JSON
dfcomb	Phase I/II Adaptive Dose-Finding Design for Combination Studies
dfcrm	Dose-finding by the continual reassessment method
dfexplore	Explore data.frames by plotting NA and classes of each variable
DFIT	Differential Functioning of Items and Tests
dfmta	Phase I/II Adaptive Dose-Finding Design for MTA
dfoptim	Derivative-free Optimization
dga	Capture-Recapture Estimation using Bayesian Model Averaging
dglars	Differential Geometric LARS (dgLARS) method
dglm	Double Generalized Linear Models
dgmb	Simulating Data for PLS Mode B Structural Models
dgof	Discrete Goodness-of-Fit Tests
dhglm	Double Hierarchical Generalized Linear Models
dHSIC	Independence Testing via Hilbert Schmidt Independence Criterion
diagonals	Block Diagonal Extraction or Replacement
diagram	Functions for visualising simple graphs (networks), plotting flow diagrams
DiagrammeR	Create Graph Diagrams and Flowcharts Using R
DiagrammeRsvg	Export DiagrammeR Graphviz Graphs as SVG
DiagTest3Grp	Diagnostic test summary measures for three ordinal groups
diaplt	Beads Summary Plot of Ranges
dice	Calculate probabilities of various dice-rolling events
DiceDesign	Designs of Computer Experiments
DiceEval	Construction and Evaluation of Metamodels
DiceKriging	Kriging Methods for Computer Experiments
DiceOptim	Kriging-Based Optimization for Computer Experiments
DiceView	Plot methods for computer experiments design and surrogate
dichromat	Color Schemes for Dichromats
dicionariosIBGE	Dictionaries for reading microdata surveys from IBGE
dielectric	Defines some physical constants and dielectric functions commonly used in optics, plasmonics
diezeit	R Interface to the ZEIT ONLINE Content API
DIFboost	Detection of Differential Item Functioning (DIF) in Rasch Models by Boosting Techniques
DiffCorr	Analyzing and Visualizing Differential Correlation Networks in Biological Data
diffdepprop	Calculates Confidence Intervals for two Dependent Proportions
diffEq	Functions from the book Solving Differential Equations in R
diffeR	Metrics of Difference for Comparing Pairs of Maps
diffIRT	Diffusion IRT Models for Response and Response Time Data
diffr	Display Differences Between Two Files using Codediff Library
diffractometry	Baseline identification and peak decomposition for x-ray diffractograms
diffusionMap	Diffusion map
DiffusionRgqd	Inference and Analysis for Generalized Quadratic Diffusions
DiffusionRjgqd	Inference and Analysis for Jump Generalized Quadratic Diffusions
DIFlasso	A penalty approach to Differential Item Functioning in Rasch Models
difR	Collection of Methods to Detect Dichotomous Differential Item Functioning (DIF)
DIFtree	Item Focused Trees for the Identification of Items in Differential Item Functioning
digest	Create Compact Hash Digests of R Objects
Digiroo2	An application programming interface for generating null models of social contacts based on individuals' space use
digitalPCR	Estimate Copy Number for Digital PCR
dils	Data-Informed Link Strength. Combine multiple-relationship networks into a single weighted network. Impute (fill-in) missing network links
DIME	DIME (Differential Identification using Mixture Ensemble)
dina	Bayesian Estimation of DINA Model
dinamic	DiNAMIC A Method To Analyze Recurrent DNA Copy Number Aberrations in Tumors
diptest	Hartigan's Dip Test Statistic for Unimodality - Corrected
DIRECT	Bayesian Clustering of Multivariate Data Under the Dirichlet-Process Prior
Directional	Directional Statistics
directlabels	Direct Labels for Multicolor Plots
directPA	Direction Analysis for Pathways and Kinases
DirichletReg	Dirichlet Regression in R
dirmult	Estimation in Dirichlet-Multinomial distribution
Disake	Discrete associated kernel estimators
disclap	Discrete Laplace Exponential Family
disclapmix	Discrete Laplace Mixture Inference using the EM Algorithm
DiscML	DiscML: An R package for estimating evolutionary rates of discrete characters using maximum likelihood
DiscreteInverseWeibull	Discrete inverse Weibull distribution
DiscreteLaplace	Discrete Laplace Distributions
discreteMTP	Multiple testing procedures for discrete test statistics
discreteRV	Create and Manipulate Discrete Random Variables
DiscreteWeibull	Discrete Weibull Distributions (Type 1 and 3)
discretization	Data preprocessing, discretization for classification
discrimARTs	Discrimination of Alternative Reproductive Tactics (ARTs)
DiscriMiner	Tools of the Trade for Discriminant Analysis
discSurv	Discrete Time Survival Analysis
diseasemapping	Modelling Spatial Variation in Disease Risk for Areal Data
DisimForMixed	Calculate Dissimilarity Matrix for Dataset with Mixed Attributes
dismo	Species Distribution Modeling
disp2D	2D Hausdorff and Simplex Dispersion Orderings
disparityfilter	Disparity Filter Algorithm of Weighted Network
displayHTS	displayHTS
dispmod	Dispersion models
disposables	Create Disposable R Packages for Testing
dissUtils	Utilities for making pairwise comparisons of multivariate data
Distance	Distance Sampling Detection Function and Abundance Estimation
distance.sample.size	Calculates Study Size Required for Distance Sampling
DistatisR	DiSTATIS Three Way Metric Multidimensional Scaling
distcomp	Computations over Distributed Data without Aggregation
distfree.cr	Distribution-free confidence region (distfree.cr)
distillery	Method Functions for Confidence Intervals and to Distill Information from an Object
distory	Distance Between Phylogenetic Histories
distr	Object oriented implementation of distributions
distrDoc	Documentation for Packages distr, distrEx, distrSim, distrTEst, distrTeach, distrMod, and distrEllipse
distrEllipse	S4 classes for elliptically contoured distributions
distrEx	Extensions of package distr
DistributionUtils	Distribution Utilities
distrMod	Object Oriented Implementation of Probability Models
distrom	Distributed Multinomial Regression
distrRmetrics	Package distr classes for distributions from Rmetrics
distrSim	Simulation classes based on package distr
distrTeach	Extensions of package distr for teaching Stochastics/Statistics in secondary school
distrTEst	Estimation and Testing classes based on package distr
divagis	Provides tools for quality checks of georeferenced plant species accessions
DivE	Diversity Estimator
diveMove	Dive Analysis and Calibration
diverse	Diversity Measures for Complex Systems
diversitree	Comparative Phylogenetic Analyses of Diversification
diveRsity	A Comprehensive, General Purpose Population Genetics Analysis Package
DiversityOccupancy	Building Diversity Models from Multiple Species Occupancy Models
DiversitySampler	Functions for re-sampling a community matrix to compute diversity indices at different sampling levels
DivMelt	HRM Diversity Assay Analysis Tool
divo	Tools for Analysis of Diversity and Similarity in Biological Systems
dixon	Nearest Neighbour Contingency Table Analysis
DJL	Distance Measure Based Judgment and Learning
dkDNA	Diffusion Kernels on a Set of Genotypes
dlm	Bayesian and Likelihood Analysis of Dynamic Linear Models
dlmap	Detection Localization Mapping for QTL
dlmodeler	Generalized Dynamic Linear Modeler
DLMtool	Data-Limited Methods Toolkit
dlnm	Distributed Lag Non-Linear Models
dma	Dynamic Model Averaging
dml	Distance Metric Learning in R
dmm	Dyadic Mixed Model for Pedigree Data
dMod	Dynamic Modeling and Parameter Estimation in ODE Models
DMR	Delete or Merge Regressors for linear model selection
dmt	Dependency Modeling Toolkit
DMwR	Functions and data for "Data Mining with R"
dna	Differential Network Analysis
DNAprofiles	DNA Profiling Evidence Analysis
DNAseqtest	Generating and Testing DNA Sequences
DNAtools	Tools for analysing forensic genetic DNA data
DnE	Distribution and Equation
dnet	Integrative Analysis of Omics Data in Terms of Network, Evolution and Ontology
DNMF	Discriminant Non-Negative Matrix Factorization
DOBAD	Analysis of Discretely Observed Linear Birth-and-Death(-and-Immigration) Markov Chains
doBy	Groupwise Statistics, LSmeans, Linear Contrasts, Utilities
docopt	Command-Line Interface Specification Language
docopulae	Optimal Designs for Copula Models
documair	Automatic Documentation for R packages
docxtractr	Extract Data Tables from Microsoft Word Documents
Dodge	Functions for Acceptance Sampling Ideas originated by H.F. Dodge
DoE.base	Full Factorials, Orthogonal Arrays and Base Utilities for DoE Packages
DoE.wrapper	Wrapper package for design of experiments functionality
doMC	Foreach Parallel Adaptor for 'parallel'
Dominance	ADI (average dominance index), social network graphs with dual directions, and music notation graph
domino	R Console Bindings for the 'Domino Command-Line Client'
doMPI	Foreach parallel adaptor for the Rmpi package
doParallel	Foreach Parallel Adaptor for the 'parallel' Package
doRedis	Foreach parallel adapter for the rredis package
doRNG	Generic Reproducible Parallel Backend for foreach Loops
DoseFinding	Planning and Analyzing Dose Finding Experiments
doSNOW	Foreach Parallel Adaptor for the 'snow' Package
dosresmeta	Performing Multivariate Dose-Response Meta-Analysis
dostats	Compute Statistics Helper Functions
dotenv	Load environment variables from .env
dotwhisker	Dot-and-Whisker Plots of Regression Results
DoubleCone	Test against parametric regression function
DoubleExpSeq	Differential Exon Usage Test for RNA-Seq Data via Empirical Bayes Shrinkage of the Dispersion Parameter
DOvalidation	Local Linear Hazard Estimation with Do-Validated and Cross-Validated Bandwidths
Dowd	Functions Ported from 'MMR2' Toolbox Offered in Kevin Dowd's Book Measuring Market Risk
downloader	Download Files over HTTP and HTTPS
downscale	Downscaling Species Occupancy
dpa	Dynamic Path Approach
dpcR	Digital PCR Analysis
dpglasso	Primal Graphical Lasso
dplR	Dendrochronology Program Library in R
dplRCon	Concordance for Dendroclimatology
dplyr	A Grammar of Data Manipulation
dpmixsim	Dirichlet Process Mixture model simulation for clustering and image segmentation
dpmr	Data Package Manager for R
DPpackage	Bayesian nonparametric modeling in R
dprep	Data Pre-Processing and Visualization Functions for Classification
dr	Methods for Dimension Reduction for Regression
drat	Drat R Archive Template
drawExpression	Visualising R syntax through graphics
DRaWR	Discriminative Random Walk with Restart
drc	Analysis of Dose-Response Curves
drfit	Dose-Response Data Evaluation
drgee	Doubly Robust Generalized Estimating Equations
DRIP	Discontinuous Regression and Image Processing
drLumi	Multiplex Immunoassays Data Analysis
drm	Regression and association models for repeated categorical data
drmdel	Dual Empirical Likelihood Inference under Density Ratio Models in the Presence of Multiple Samples
dropR	Analyze Drop Out of an Experiment or Survey
drsmooth	Dose-Response Modeling with Smoothing Splines
DrugClust	Implementation of a Machine Learning Framework for Predicting Drugs Side Effects
ds	Descriptive Statistics
dsample	Discretization-Based Direct Random Sample Generation
DSBayes	Bayesian subgroup analysis in clinical trials
dse	Dynamic Systems Estimation (Time Series Package)
DSL	Distributed Storage and List
dslice	Dynamic Slicing
dsm	Density Surface Modelling of Distance Sampling Data
DSpat	Spatial Modelling for Distance Sampling Data
DSsim	Distance Sampling Simulations
dst	Using Dempster-Shafer Theory
DStree	Recursive Partitioning for Discrete-Time Survival Trees
DSviaDRM	Exploring Disease Similarity in Terms of Dysfunctional Regulatory Mechanisms
DT	A Wrapper of the JavaScript Library 'DataTables'
DTComPair	Comparison of Binary Diagnostic Tests in a Paired Study Design
DTDA	Doubly truncated data analysis
dti	Analysis of Diffusion Weighted Imaging (DWI) Data
DTK	Dunnett-Tukey-Kramer Pairwise Multiple Comparison Test Adjusted for Unequal Variances and Unequal Sample Sizes
DTMCPack	Suite of functions related to discrete-time discrete-state Markov Chains
DTR	Estimation and Comparison of Dynamic Treatment Regimes
DTRlearn	Learning Algorithms for Dynamic Treatment Regimes
DTRreg	DTR Estimation and Inference via G-Estimation, Dynamic WOLS, and Q-Learning
dtt	Discrete Trigonometric Transforms
dtw	Dynamic Time Warping Algorithms
dtwclust	Time Series Clustering Along with Optimizations for the Dynamic Time Warping Distance
dtwSat	Time-Weighted Dynamic Time Warping for Remote Sensing Time Series Analysis
dualScale	Dual Scaling Analysis of Multiple Choice Data
dummies	Create dummy/indicator variables flexibly and efficiently
dummy	Automatic Creation of Dummies with Support for Predictive Modeling
DunnettTests	Software implementation of step-down and step-up Dunnett test procedures
dunn.test	Dunn's Test of Multiple Comparisons Using Rank Sums
dupiR	Bayesian inference from count data using discrete uniform priors
dvfBm	Discrete variations of a fractional Brownian motion
DVHmetrics	Analyze Dose-Volume Histograms and Check Constraints
dvn	Access to The Dataverse Network APIs
DWreg	Parametric Regression for Discrete Response
dygraphs	Interface to 'Dygraphs' Interactive Time Series Charting Library
DYM	Did You Mean?
dyn	Time Series Regression
DynamicDistribution	Dynamically visualized probability distributions and their moments
dynamicGraph	dynamicGraph
dynamicTreeCut	Methods for Detection of Clusters in Hierarchical Clustering Dendrograms
dynatopmodel	Implementation of the Dynamic TOPMODEL Hydrological Model
dynaTree	Dynamic trees for learning and design
dynBiplotGUI	Full Interactive GUI for Dynamic Biplot in R
DynClust	Denoising and clustering for dynamical image sequence (2D or 3D)+T
dynCorr	Dynamic Correlation Package
dynia	Fit Dynamic Intervention Model
dynlm	Dynamic Linear Regression
DynNom	Dynamic Nomograms for Linear, Generalized Linear and Proportional Hazard Models
dynpred	Companion Package to "Dynamic Prediction in Clinical Survival Analysis"
dynRB	Dynamic Range Boxes
dynsim	Dynamic Simulations of Autoregressive Relationships
dynsurv	Dynamic models for survival data
DynTxRegime	Methods for Estimating Dynamic Treatment Regimes
e1071	Misc Functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien
eaf	Plots of the Empirical Attainment Function
earlywarnings	Early Warning Signals Toolbox for Detecting Critical Transitions in Timeseries
earth	Multivariate Adaptive Regression Splines
EasyABC	Efficient Approximate Bayesian Computation Sampling Schemes
easyanova	Analysis of variance and other important complementary analyzes
EasyHTMLReport	EasyHTMLReport
EasyMARK	Utility functions for working with mark-recapture data
easynls	Easy nonlinear model
easypackages	Easy Loading and Installing of Packages
easypower	Sample Size Estimation for Experimental Designs
easyPubMed	Search and Retrieve Scientific Publication Records from Pubmed
EasyStrata	Evaluation of stratified genome-wide association meta-analysis results
easyVerification	Ensemble Forecast Verification for Large Data Sets
eba	Elimination-by-Aspects (EBA) Models
ebal	Entropy reweighting to create balanced samples
EbayesThresh	Empirical Bayes Thresholding and Related Methods
ebdbNet	Empirical Bayes Estimation of Dynamic Bayesian Networks
EBEN	Empirical Bayesian Elastic Net
ebGenotyping	Genotyping using Next Generation Sequencing Data
EBglmnet	Empirical Bayesian Lasso and Elastic Net Methods for Generalized Linear Models
EBMAforecast	Ensemble BMA Forecasting
EBS	Exact Bayesian Segmentation
ebSNP	Genotyping and SNP calling using single-sample next generation sequencing data
ecb	Programmatic Access to the European Central Bank's Statistical Data Warehouse (SDW)
ECctmc	Simulation from Endpoint-Conditioned Continuous Time Markov Chains
ecd	Elliptic Distribution Based on Elliptic Curves
Ecdat	Data Sets for Econometrics
ecespa	Functions for Spatial Point Pattern Analysis
Ecfun	Functions for Ecdat
ecipex	Efficient calculation of fine structure isotope patterns via Fourier transforms of simplex-based elemental models
eco	Ecological Inference in 2x2 Tables
ecodist	Dissimilarity-based functions for ecological analysis
ecoengine	Programmatic Interface to the API Serving UC Berkeley's Natural History Data
EcoGenetics	Spatial Analysis of Phenotypic, Genotypic and Environmental Data
EcoHydRology	A community modeling foundation for Eco-Hydrology
ecolMod	"A practical guide to ecological modelling - using R as a simulation platform"
ecoreg	Ecological Regression using Aggregate and Individual Data
ecoretriever	R Interface to the EcoData Retriever
ecosim	Toolbox for Aquatic Ecosystem Modeling
EcoSimR	Null Model Analysis for Ecological Data
ECOSolveR	Embedded Conic Solver in R
ecospace	Simulating Community Assembly and Ecological Diversification Using Ecospace Frameworks
ecospat	Spatial Ecology Miscellaneous Methods
ecotoxicology	Methods for Ecotoxicology
EcoTroph	EcoTroph R package
ecoval	Procedures for Ecological Assessment of Surface Waters
EcoVirtual	Simulation of Ecological Models
ecp	Non-Parametric Multiple Change-Point Analysis of Multivariate Data
edcc	Economic Design of Control Charts
edeaR	Exploratory and Descriptive Event-Based Data Analysis
edeR	Email Data Extraction Using R
edesign	Maximum Entropy Sampling
EDFIR	Estimating Discrimination Factors
edfReader	Reading EDF(+) and BDF(+) Files
edgar	Platform for EDGAR Filing Management
edgebundleR	Circle Plot with Bundled Edges
edgeCorr	Spatial Edge Correction
edgeRun	More Powerful Unconditional Testing of Negative Binomial Means for Digital Gene Expression Data
EDISON	Network Reconstruction and Changepoint Detection
EditImputeCont	Simultaneous Edit-Imputation for Continuous Microdata
editrules	Parsing, Applying, and Manipulating Data Cleaning Rules
edmr	Empirical Differentially Methylated Regions Calculation
EDR	Estimation of the effective dimension reduction (EDR) space
edrGraphicalTools	Provides tools for dimension reduction methods
eegAnalysis	Tools for analysis and classification of electroencephalography (EEG) data
eegkit	Toolkit for Electroencephalography Data
eegkitdata	Data for package eegkit
eel	Extended Empirical Likelihood
EEM	Read and Preprocess Fluorescence Excitation-Emission Matrix (EEM) Data
eemR	Tools for Pre-Processing Emission-Excitation-Matrix (EEM) Fluorescence Data
eeptools	Convenience Functions for Education Data
EFDR	Wavelet-Based Enhanced FDR for Signal Detection in Noisy Images
EffectLiteR	Average and Conditional Effects
effects	Effect Displays for Linear, Generalized Linear, and Other Models
EffectsRelBaseline	Test changes of a grouped response relative to baseline
EffectStars	Visualization of Categorical Response Models
EffectTreat	Prediction of Therapeutic Success
efflog	The Causal Effects for a Causal Loglinear Model
effsize	Efficient Effect Size Computation
efreadr	Read European Fluxes CSV Files
ega	Error Grid Analysis
egcm	Engle-Granger Cointegration Models
eggCounts	Hierarchical Modelling of Faecal Egg Counts
egonet	Tool for ego-centric measures in Social Network Analysis
EGRET	Exploration and Graphics for RivEr Trends (EGRET)
EGRETci	Exploration and Graphics for RivEr Trends (EGRET) Confidence Intervals
eha	Event History Analysis
eHOF	Extended HOF (Huisman-Olff-Fresco) Models
ei	Ecological Inference
EIAdata	R Wrapper for the Energy Information Administration (EIA) API
eiCompare	Compares EI, Goodman, RxC Estimates
eigeninv	Generates (dense) matrices that have a given set of eigenvalues
eigenmodel	Semiparametric factor and regression models for symmetric relational data
eigenprcomp	Computes confidence intervals for principal components
EILA	Efficient Inference of Local Ancestry
eiPack	eiPack: Ecological Inference and Higher-Dimension Data Management
eive	An algorithm for reducing errors-in-variable bias in simple linear regression
eiwild	Ecological Inference with individual and aggregate data
EL	Two-sample Empirical Likelihood
elasso	Enhanced Least Absolute Shrinkage and Selection Operator Regression Model
elastic	General Purpose Interface to 'Elasticsearch'
elasticnet	Elastic-Net for Sparse Estimation and Sparse PCA
elec	Collection of functions for statistical election audits
elec.strat	Functions for election audits using stratified random samples
ElemStatLearn	Data Sets, Functions and Examples from the Book: "The Elements of Statistical Learning, Data Mining, Inference, and Prediction" by Trevor Hastie, Robert Tibshirani and Jerome Friedman
elexr	Load Associated Press Election Results with Elex
elliplot	Ellipse Summary Plot of Quantiles
ellipse	Functions for drawing ellipses and ellipse-like confidence regions
elliptic	elliptic functions
elmNN	Implementation of ELM (Extreme Learning Machine ) algorithm for SLFN ( Single Hidden Layer Feedforward Neural Networks )
ELMR	Extreme Machine Learning (ELM)
EloChoice	Preference Rating for Visual Stimuli Based on Elo Ratings
EloRating	Animal Dominance Hierarchies by Elo Rating
elrm	Exact Logistic Regression via MCMC
ElstonStewart	Elston-Stewart Algorithm
ELT	Experience Life Tables
ELYP	Empirical Likelihood Analysis for the Cox Model and Yang-Prentice (2005) Model
EMA	Easy Microarray data Analysis
EMbC	Expectation-Maximization Binary Clustering
embryogrowth	Tools to Analyze the Thermal Reaction Norm of Embryo Growth
EMC	Evolutionary Monte Carlo (EMC) algorithm
EMCluster	EM Algorithm for Model-Based Clustering of Finite Mixture Gaussian Distribution
EMD	Empirical Mode Decomposition and Hilbert Spectral Analysis
emdbook	Support Functions and Data for "Ecological Models and Data"
emdist	Earth Mover's Distance
emg	Exponentially Modified Gaussian (EMG) Distribution
emil	Evaluation of Modeling without Information Leakage
emIRT	EM Algorithms for Estimating Item Response Theory Models
emma	Evolutionary model-based multiresponse approach
EMMAgeo	End-Member Modelling of Grain-Size Data
emme2	Read and Write to an EMME/2 databank
EMMIXcontrasts	Contrasts in mixed effects for EMMIX model with random effects
EMMIXskew	The EM Algorithm and Skew Mixture Distribution
EMMIXuskew	Fitting Unrestricted Multivariate Skew t Mixture Models
EMMREML	Fitting Mixed Models with Known Covariance Structures
emoa	Evolutionary Multiobjective Optimization Algorithms
emojifont	Emoji Fonts for using in R
emov	Eye Movement Analysis Package for Fixation and Saccade Detection
EMP	Expected Maximum Profit Classification Performance Measure
EmpiricalCalibration	Routines for Performing Empirical Calibration of Observational Study Estimates
empiricalFDR.DESeq2	Simulation-Based False Discovery Rate in RNA-Seq
emplik	Empirical Likelihood Ratio for Censored/Truncated Data
emplik2	Empirical Likelihood Ratio Test for Two Samples with Censored Data
EMT	Exact Multinomial Test: Goodness-of-Fit Test for Discrete Multivariate data
emulator	Bayesian emulation of computer programs
emuR	Main Package of the EMU Speech Database Management System
EMVC	Entropy Minimization over Variable Clusters (EMVC)
enaR	Tools for Ecological Network Analysis
endogMNP	R Package for Fitting Multinomial Probit Models with Endogenous Selection
endorse	R Package for Analyzing Endorsement Experiments
energy	E-statistics (energy statistics)
english	Translate integers into English
EngrExpt	Data sets from "Introductory Statistics for Engineering Experimentation"
enigma	Client for the 'Enigma' 'API'
ENiRG	Ecological Niche in R and GRASS
ENMeval	Automated Runs and Evaluations of Ecological Niche Models
ENmisc	Neuwirth miscellaneous
enpls	Ensemble Partial Least Squares Regression
EnQuireR	A package dedicated to questionnaires
enRich	An R Package for the Analysis of Multiple ChIP-Seq Data
enrichvs	Enrichment assessment of virtual screening approaches
EnsembleBase	Extensible Package for Parallel, Batch Training of Base Learners for Ensemble Modeling
ensembleBMA	Probabilistic Forecasting using Ensembles and Bayesian Model Averaging
EnsembleCV	Extensible Package for Cross-Validation-Based Integration of Base Learners
ensembleMOS	Ensemble Model Output Statistics
EnsemblePCReg	Extensible Package for Principal-Component-Regression-Based Integration of Base Learners
EnsemblePenReg	Extensible Classes and Methods for Penalized-Regression-based Integration of Base Learners
ensurer	Ensure Values at Runtime
entropart	Entropy Partitioning to Measure Diversity
entropy	Estimation of Entropy, Mutual Information and Related Quantities
EntropyEstimation	Estimation of Entropy and Related Quantities
EntropyExplorer	Tools for Exploring Differential Shannon Entropy, Differential Coefficient of Variation and Differential Expression
enveomics.R	Various R Functions from the Kostas Lab
enviPat	Isotope Pattern, Profile and Centroid Calculation for Mass Spectrometry
enviPick	Peak Picking for High Resolution Mass Spectrometry Data
EnviroStat	Statistical Analysis of Environmental Space-Time Processes
envlpaster	Enveloping the Aster Model
EnvNicheR	Niche Estimation
EnvStats	Package for Environmental Statistics, Including US EPA Guidance
epade	Easy Plots
epandist	Statistical Functions for the Censored and Uncensored Epanechnikov Distribution
epanetReader	Read Epanet Files into R
EPGLM	Gaussian Approximation of Bayesian Binary Regression Models
Epi	A Package for Statistical Analysis in Epidemiology
epibasix	Elementary Epidemiological Functions for Epidemiology and Biostatistics
EpiBayes	Implements Hierarchical Bayesian Models for Epidemiological Applications
EpiContactTrace	Epidemiological Tool for Contact Tracing
epiDisplay	Epidemiological Data Display Package
EpiDynamics	Dynamic Models in Epidemiology
EpiEstim	EpiEstim: a package to estimate time varying reproduction numbers from epidemic curves
epifit	Flexible Modelling Functions for Epidemiological Data Analysis
EpiModel	Mathematical Modeling of Infectious Disease
epinet	Epidemic/Network-Related Tools
epiR	Tools for the Analysis of Epidemiological Data
episensr	Basic Sensitivity Analysis of Epidemiological Results
episplineDensity	Density Estimation with Soft Information by Exponential Epi-splines
epitools	Epidemiology Tools
Eplot	Plotting longitudinal series
epoc	EPoC (Endogenous Perturbation analysis of Cancer)
epr	Easy polynomial regression
EQL	Extended-Quasi-Likelihood-Function (EQL)
eqs2lavaan	EQS Output Conversion to lavaan Functions
eqtl	Tools for analyzing eQTL experiments: A complementary to Karl Broman's 'qtl' package for genome-wide analysis
equate	Observed-Score Linking and Equating
equateIRT	Direct, Chain and Average Equating Coefficients with Standard Errors Using IRT Methods
equivalence	Provides Tests and Graphics for Assessing Tests of Equivalence
erboost	Nonparametric Multiple Expectile Regression via ER-Boost
erer	Empirical Research in Economics with R
ergm	Fit, Simulate and Diagnose Exponential-Family Models for Networks
ergm.count	Fit, Simulate and Diagnose Exponential-Family Models for Networks with Count Edges
ergm.ego	Fit, Simulate and Diagnose Exponential-Family Random Graph Models to Egocentrically Sampled Network Data
ergm.graphlets	ERG Modeling Based on Graphlet Properties
ergmharris	Local Health Department network data set
ergm.userterms	User-specified terms for the statnet suite of packages
eRm	Extended Rasch Modeling
ERP	Significance Analysis of Event-Related Potentials Data
erp.easy	Event-Related Potential (ERP) Data Exploration Made Easy
erpR	Event-related potentials (ERP) analysis, graphics and utility functions
errint	Build Error Intervals
ES	Edge Selection
esaBcv	Estimate Number of Latent Factors and Factor Matrix for Factor Analysis
ESEA	ESEA: Discovering the Dysregulated Pathways based on Edge Set Enrichment Analysis
ESG	ESG - A package for asset projection
ESGtoolkit	Toolkit for the simulation of financial assets and interest rates models
ESKNN	Ensemble of Subset of K-Nearest Neighbours Classifiers for Classification and Class Membership Probability Estimation
estatapi	R Interface to e-Stat API
EstCRM	Calibrating Parameters for the Samejima's Continuous IRT Model
EstHer	Estimation of Heritability in High Dimensional Sparse Linear Mixed Models using Variable Selection
estimability	Tools for Assessing Estimability of Linear Predictions
estout	Estimates Output
EstSimPDMP	Estimation and Simulation for PDMPs
etable	Easy Table
ETAS	Modeling Earthquake Data Using ETAS Model
etasFLP	Mixed FLP and ML Estimation of ETAS Space-Time Point Processes
ETC	Equivalence to control
ETLUtils	Utility Functions to Execute Standard Extract/Transform/Load Operations (using Package 'ff') on Large Data
etm	Empirical Transition Matrix
etma	Epistasis Test in Meta-Analysis
eulerian	eulerian: A package to find eulerian paths from graphs
euroMix	Calculations for DNA Mixtures
EurosarcBayes	Bayesian Single Arm Sample Size Calculation Software
eurostat	Tools for Eurostat Open Data
eva	Extreme Value Analysis with Goodness-of-Fit Testing
EvalEst	Dynamic Systems Estimation - Extensions
evaluate	Parsing and Evaluation Tools that Provide More Details than the Default
Evapotranspiration	Modelling Actual, Potential and Reference Crop Evapotranspiration
EvCombR	Evidence Combination in R
evd	Functions for Extreme Value Distributions
evdbayes	Bayesian Analysis in Extreme Value Theory
eVenn	A Powerful Tool to Quickly Compare Huge Lists and Draw Venn Diagrams
eventInterval	Sequential Event Interval Analysis
events	Store and manipulate event data
eventstudies	Event study and extreme event analysis
evir	Extreme Values in R
evmix	Extreme Value Mixture Modelling, Threshold Estimation and Boundary Corrected Kernel Density Estimation
evobiR	Comparative and Population Genetic Analyses
evolqg	Tools for Evolutionary Quantitative Genetics
evolvability	Calculation of Evolvability Parameters
Evomorph	Evolutionary Morphometric Simulation
EvoRAG	Evolutionary Rates Across Gradients
evt0	Mean of order p, peaks over random threshold Hill and high quantile estimates
evtree	Evolutionary Learning of Globally Optimal Trees
EW	Edgeworth Expansion
EWGoF	Goodness-of-Fit Tests for the Exponential and Two-Parameter Weibull Distributions
Exact	Unconditional Exact Test
exact2x2	Exact Conditional Tests and Confidence Intervals for 2x2 tables
exactci	Exact P-Values and Matching Confidence Intervals for Simple Discrete Parametric Cases
ExactCIdiff	Inductive Confidence Intervals for the difference between two proportions
exactLoglinTest	Monte Carlo Exact Tests for Log-linear models
exactmeta	Exact fixed effect meta analysis
ExactPath	Exact solution paths for regularized LASSO regressions with L_1 penalty
exactRankTests	Exact Distributions for Rank and Permutation Tests
exams	Automatic Generation of Exams in R
ExceedanceTools	Confidence regions for exceedance sets and contour lines
excel.link	Convenient Data Exchange with Microsoft Excel
exCon	Interactive Exploration of Contour Data
excursions	Excursion Sets and Contour Credibility Regions for Random Fields
exif	Read EXIF Metadata from JPEGs
exifr	EXIF Image Data in R
ExomeDepth	Calls Copy Number Variants from Targeted Sequence Data
expands	Expanding Ploidy and Allele-Frequency on Nested Subpopulations
ExpDE	Modular Differential Evolution for Experimenting with Operators
ExpDes	Experimental Designs package
ExpDes.pt	Pacote Experimental Designs (Portuguese)
expectreg	Expectile and Quantile Regression
experiment	experiment: R package for designing and analyzing randomized experiments
expert	Modeling without data using expert opinion
ExplainPrediction	Explanation of Predictions for Classification and Regression Models
explor	Interactive Interfaces for Results Exploration
exploreR	Tools for Quickly Exploring Data
expm	Matrix Exponential
expoRkit	Expokit in R
ExPosition	Exploratory analysis with the singular value decomposition
expoTree	Calculate density dependent likelihood of a phylogenetic tree
expp	Spatial analysis of extra-pair paternity
expsmooth	Data Sets from "Forecasting with Exponential Smoothing"
exptest	Tests for Exponentiality
exreport	Fast, Reliable and Elegant Reproducible Research
exsic	Convenience Functions for Botanists to Create Specimens Indices
ExtDist	Extending the Range of Functions for Probability Distributions
extfunnel	Additional Funnel Plot Augmentations
extlasso	Maximum penalized likelihood estimation with extended lasso penalty
extraBinomial	Extra-binomial approach for pooled sequencing data
extracat	Categorical Data Analysis and Visualization
extrafont	Tools for using fonts
extrafontdb	Package for holding the database for the extrafont package
extraTrees	Extremely Randomized Trees (ExtraTrees) Method for Classification and Regression
ExtremeBounds	Extreme Bounds Analysis (EBA)
extRemes	Extreme Value Analysis
extremevalues	Univariate Outlier Detection
extremogram	Estimation of Extreme Value Dependence for Time Series Data
extWeibQuant	Estimate Lower Extreme Quantile with the Censored Weibull MLE and Censored Weibull Mixture
eyetracking	Eyetracking Helper Functions
eyetrackingR	Eye-Tracking Data Analysis
ez	Easy Analysis and Visualization of Factorial Experiments
ezec	Easy Interface to Effective Concentration Calculations
ezglm	selects significant non-additive interaction between two variables using fast GLM implementation
ezknitr	Avoid the Typical Working Directory Pain When Using 'knitr'
ezsim	provide an easy to use framework to conduct simulation
ezsummary	Summarise Data in the Quick and Easy Way
FacPad	Bayesian Sparse Factor Analysis model for the inference of pathways responsive to drug treatment
FactMixtAnalysis	Factor Mixture Analysis with covariates
FACTMLE	Maximum Likelihood Factor Analysis
FactoClass	Combination of Factorial Methods and Cluster Analysis
FactoMineR	Multivariate Exploratory Data Analysis and Data Mining
factorplot	factorplot
factorQR	Bayesian quantile regression factor models
Factoshiny	Perform Factorial Analysis from FactoMineR with a Shiny Application
FACTscorer	Scores the FACT and FACIT Family of Patient-Reported Outcome Measures
factualR	thin wrapper for the Factual.com server API
FADA	Variable selection for supervised classification in high dimension
FAdist	Distributions that are Sometimes Used in Hydrology
Fahrmeir	Data from the Book "Multivariate Statistical Modelling Based on Generalized Linear Models", First Edition, by Ludwig Fahrmeir and Gerhard Tutz
fail	File Abstraction Interface Layer (FAIL)
FAiR	Factor Analysis in R
faisalconjoint	Faisal Conjoint Model: A New Approach to Conjoint Analysis
falcon	Finding Allele-specific Copy Number in Next-Generation Sequencing Data
falsy	Define Truthy and Falsy Values
fame	Interface for FAME Time Series Database
FamEvent	Simulation of Time-to-Event Family Data and Penetrance Estimation
Familias	Probabilities for Pedigrees Given DNA Data
FAMILY	A Convex Formulation for Modeling Interactions with Strong Heredity
FAmle	Maximum Likelihood and Bayesian Estimation of Univariate Probability Distributions
FAMT	Factor Analysis for Multiple Testing (FAMT) : simultaneous tests under dependence in high-dimensional data
fanc	Penalized Likelihood Factor Analysis via Nonconvex Penalty
fANCOVA	Nonparametric Analysis of Covariance
fancycut	A Fancy Version of 'base::cut'
fanovaGraph	Building Kriging Models from FANOVA Graphs
fanplot	Visualisation of Sequential Probability Distributions Using Fan Charts
FAOSTAT	Download Data from the FAOSTAT Database of the Food and Agricultural Organization (FAO) of the United Nations
faoutlier	Influential Case Detection Methods for Factor Analysis and Structural Equation Models
far	Modelization for Functional AutoRegressive Processes
faraway	Functions and Datasets for Books by Julian Faraway
fArma	ARMA Time Series Modelling
farsi	Translate integers into persian
fAsianOptions	EBM and Asian Option Valuation
fAssets	Rmetrics - Analysing and Modelling Financial Assets
fast	Implementation of the Fourier Amplitude Sensitivity Test (FAST)
fastAdaboost	a Fast Implementation of Adaboost
FastBandChol	Fast Estimation of a Covariance Matrix by Banding the Cholesky Factor
fastclime	A Fast Solver for Parameterized Linear Programming Problems and Constrained L1 Minimization Approach to Sparse Precision Matrix Estimation
fastcluster	Fast Hierarchical Clustering Routines for R and Python
fastcox	Lasso and elastic-net penalized Cox's regression in high dimensions models using the cocktail algorithm
fastdigest	Fast, Low Memory-Footprint Digests of R Objects
fastGHQuad	Fast Rcpp implementation of Gauss-Hermite quadrature
FastGP	Efficiently Using Gaussian Processes with Rcpp and RcppEigen
fastGraph	Fast Drawing and Shading of Graphs of Statistical Distributions
FastHCS	Robust Algorithm for Principal Component Analysis
fastHICA	Hierarchical Independent Component Analysis: a Multi-Scale Sparse Non-Orthogonal Data-Driven Basis
fastICA	FastICA Algorithms to perform ICA and Projection Pursuit
FastImputation	Learn from Training Data then Quickly Fill in Missing Data
FastKM	A Fast Multiple-Kernel Method Based on a Low-Rank Approximation
FastKNN	Fast k-Nearest Neighbors
fastM	Fast Computation of Multivariate M-estimators
fastmatch	Fast match() function
FastPCS	FastPCS Robust Fit of Multivariate Location and Scatter
fastpseudo	Fast Pseudo Observations
fastR	Foundations and Applications of Statistics Using R
FastRCS	Fits the FastRCS Robust Multivariable Linear Regression Model
FastRWeb	Fast Interactive Framework for Web Scripting Using R
fastSOM	Fast Calculation of Spillover Measures
fasttime	Fast Utility Function for Time Parsing and Conversion
fat2Lpoly	Two-Locus Family-Based Association Test with Polytomic Outcome
FatTailsR	Kiener Distributions and Fat Tails in Finance
favnums	A Dataset of Favourite Numbers
FAwR	Functions and Datasets for "Forest Analytics with R"
fBasics	Rmetrics - Markets and Basic Statistics
fbati	Gene by Environment Interaction and Conditional Gene Tests for Nuclear Families
FBFsearch	Algorithm for searching the space of Gaussian directed acyclic graphical models through moment fractional Bayes factors
FBN	FISH Based Normalization and Copy Number inference of SNP microarray data
fBonds	Bonds and Interest Rate Models
fbRanks	Association Football (Soccer) Ranking via Poisson Regression
fbroc	Fast Algorithms to Bootstrap Receiver Operating Characteristics Curves
fcd	Fused Community Detection
fCertificates	Basics of Certificates and Structured Products Valuation
FCGR	Fatigue Crack Growth in Reliability
fclust	Fuzzy Clustering
FCMapper	Fuzzy Cognitive Mapping
FCNN4R	Fast Compressed Neural Networks for R
fCopulae	Rmetrics - Bivariate Dependence Structures with Copulae
fcros	A Method to Search for Differentially Expressed Genes and to Detect Recurrent Chromosomal Copy Number Aberrations
FD	Measuring functional diversity (FD) from multiple traits, and other tools for functional ecology
fda	Functional Data Analysis
fdakma	Functional Data Analysis: K-Mean Alignment
fdaMixed	Functional data analysis in a mixed model framework
fdapace	Functional Data Analysis and Empirical Dynamics
fdaPDE	Functional Data Analysis and Partial Differential Equations; Statistical Analysis of Functional and Spatial Data, Based on Regression with Partial Differential Regularizations
fdasrvf	Elastic Functional Data Analysis
fdatest	Interval Testing Procedure for Functional Data
fda.usc	Functional Data Analysis and Utilities for Statistical Computing
FDboost	Boosting Functional Regression Models
FDGcopulas	Multivariate Dependence with FDG Copulas
fdrci	Permutation-based FDR Point and Confidence Interval Estimation
fdrDiscreteNull	False Discovery Rate Procedure Under Discrete Null Distributions
FDRreg	False discovery rate regression
FDRsampsize	Compute Sample Size that Meets Requirements for Average Power and FDR
fdrtool	Estimation of (Local) False Discovery Rates and Higher Criticism
fds	Functional data sets
fdth	Frequency Distribution Tables, Histograms and Polygons
FeaLect	Scores Features for Feature Selection
feature	Local Inferential Feature Significance for Multivariate Kernel Density Estimation
FeatureHashing	Creates a Model Matrix via Feature Hashing with a Formula Interface
features	Feature Extraction for Discretely-Sampled Functional Data
fechner	Fechnerian Scaling of Discrete Object Sets
FedData	Functions to Automate Downloading Geospatial Data Available from Several Federated Data Sources
federalregister	Client Package for the U.S. Federal Register API
FeedbackTS	Analysis of Feedback in Time Series
FENmlm	Fixed Effects Nonlinear Maximum Likelihood Models
fermicatsR	Fermi Large Area Telescope Catalogs
fExoticOptions	Exotic Option Valuation
fExpressCertificates	fExpressCertificates - Structured Products Valuation for ExpressCertificates/Autocallables
fExtremes	Rmetrics - Extreme Financial Market Data
ff	memory-efficient storage of large data on disk and fast access functions
ffbase	Basic Statistical Functions for Package 'ff'
FFD	Freedom from Disease
FField	Force field simulation for a set of points
ffmanova	Fifty-fifty MANOVA
fftw	Fast FFT and DCT based on FFTW
fftwtools	Wrapper for FFTW3: Includes 1-D, Univariate and Multivariate, and 2-D Transform
fgac	Generalized Archimedean Copula
FGalgorithm	Flury and Gautschi algorithms
fGarch	Rmetrics - Autoregressive Conditional Heteroskedastic Modelling
Fgmutils	Forest Growth Model Utilities
FGN	Fractional Gaussian Noise and power law decay time series model fitting
fgof	Fast Goodness-of-fit Test
fgpt	Floating Grid Permutation Technique
FGSG	Feature Grouping and Selection Over an Undirected Graph
fgui	Function GUI
fheatmap	Fantastic Heatmap
FHtest	Tests for Right and Interval-Censored Survival Data Based on the Fleming-Harrington Class
FI	Provide functions for forest inventory calculations
FIACH	Retrospective Noise Control for fMRI
fICA	Classical, Reloaded and Adaptive FastICA Algorithms
fields	Tools for Spatial Data
FieldSim	Random Fields (and Bridges) Simulations
fifer	A collection of miscellaneous functions
filehash	Simple Key-Value Database
filehashSQLite	Simple key-value database using SQLite
filematrix	File-Backed Matrix Class with Convenient Read and Write Access
filenamer	Easy Management of File Names
fImport	Rmetrics - Economic and Financial Data Import
financial	Solving financial problems in R
FinancialInstrument	Financial Instrument Model Infrastructure for R
FinAsym	Classifies implicit trading activity from market quotes and computes the probability of informed trading
FinCal	Time Value of Money, Time Series Analysis and Computational Finance
FinCovRegularization	Covariance Matrix Estimation and Regularization for Finance
FindAllRoots	Find all root(s) of the equation and Find root(s) of the equation by dichotomy
FindIt	Finding Heterogeneous Treatment Effects
FindMinIC	Find Models with Minimum IC
findpython	Python tools to find an acceptable python binary
fingerprint	Functions to operate on binary fingerprint data
finiteruinprob	Computation of the probability of ruin within a finite time horizon
FinTS	Companion to Tsay (2005) Analysis of Financial Time Series
FisherEM	The Fisher-EM algorithm
fisheyeR	Fisheye and Hyperbolic-space-alike Interactive Visualization Tools in R
FisHiCal	Iterative FISH-based Calibration of Hi-C Data
fishmethods	Fishery Science Methods and Models in R
fishMod	Fits Poisson-sum-of-Gammas GLMs, Tweedie GLMs, and delta log-normal mdoels
fishmove	Prediction of Fish Movement Parameters
fit4NM	NONMEM platform
FitAR	Subset AR Model Fitting
FitARMA	FitARMA: Fit ARMA or ARIMA using fast MLE algorithm
fitbitScraper	Scrapes Data from Fitbit
fitdistrplus	Help to Fit of a Parametric Distribution to Non-Censored or Censored Data
fitDRC	Fitting Density Ratio Classes
fit.models	fit.models
FITSio	FITS (Flexible Image Transport System) utilities
fitTetra	fitTetra is an R package for assigning tetraploid genotype scores
FKF	Fast Kalman Filter
flacco	Feature-Based Landscape Analysis of Continuous and Constraint Optimization Problems
flam	Fits Piecewise Constant Models with Data-Adaptive Knots
flare	Family of Lasso Regression
flashClust	Implementation of optimal hierarchical clustering
flexclust	Flexible Cluster Algorithms
flexCWM	Flexible Cluster-Weighted Modeling
flexmix	Flexible Mixture Modeling
FlexParamCurve	Tools to Fit Flexible Parametric Curves
flexPM	Flexible Parametric Models for Censored and Truncated Data
flexsurv	Flexible Parametric Survival and Multi-State Models
FLIM	Farewell’s Linear Increments Model
flip	Multivariate Permutation Tests
FLLat	Fused Lasso Latent Feature Model
flora	Tools for Interacting with the Brazilian Flora 2020
flowDiv	Cytometric Diversity Indices from 'FlowJo' Workspaces
flower	Tools for characterizing flowering traits
flowfield	Forecasts future values of a univariate time series
flowr	Streamlining Design and Deployment of Complex Workflows
flows	Flow Selection and Analysis
FlowScreen	Daily Streamflow Trend and Change Point Screening
FLR	Fuzzy Logic Rule Classifier
flsa	Path algorithm for the general Fused Lasso Signal Approximator
FLSSS	Fixed Size Subset Sum Solution
Flury	Data Sets from Flury, 1997
flux	Flux rate calculation from dynamic closed chamber measurements
fma	Data sets from "Forecasting: methods and applications" by Makridakis, Wheelwright & Hyndman (1998)
FME	A Flexible Modelling Environment for Inverse Modelling, Sensitivity, Identifiability, Monte Carlo Analysis
FMP	Filtered Monotonic Polynomial IRT Models
fmri	Analysis of fMRI experiments
fmrs	Variable Selection in Finite Mixture of AFT Regression and FMR
fmsb	Functions for Medical Statistics Book with some Demographic Data
FMStable	Finite Moment Stable Distributions
fmt	Variance estimation of FMT method (Fully Moderated t-statistic)
fMultivar	Rmetrics - Analysing and Modeling Multivariate Financial Return Distributions
FNN	Fast Nearest Neighbor Search Algorithms and Applications
fNonlinear	Nonlinear and Chaotic Time Series Modelling
foba	greedy variable selection
fontcm	Computer Modern font for use with extrafont package
foodweb	visualisation and analysis of food web networks
fOptions	Rmetrics - Pricing and Evaluating Basic Options
forams	Foraminifera and Community Ecology Analyses
foreach	Provides Foreach Looping Construct for R
ForeCA	Forecastable Component Analysis
forecast	Forecasting Functions for Time Series and Linear Models
ForecastCombinations	Forecast Combinations
forecTheta	Forecasting Time Series by Theta Models
forega	Floating-Point Genetic Algorithms with Statistical Forecast Based Inheritance Operator
foreign	Read Data Stored by Minitab, S, SAS, SPSS, Stata, Systat, Weka, dBase, ...
forensic	Statistical Methods in Forensic Genetics
forensim	Statistical tools for the interpretation of forensic DNA mixtures
forestFloor	Visualizes Random Forests with Feature Contributions
forestmodel	Forest Plots from Regression Models
forestplot	Advanced Forest Plot Using 'grid' Graphics
ForImp	Imputation of Missing Values Through a Forward Imputation Algorithm
ForIT	Functions from the 2nd Italian Forest Inventory (INFC)
FormalSeries	Elementary arithemtic in formal series rings
formatR	Format R Code Automatically
formattable	Formattable Data Structures
Formula	Extended Model Formulas
formula.tools	Utilities for Formulas, Expressions, Calls and Other Objects
fortunes	R Fortunes
forward	Forward search
ForwardSearch	Forward Search using asymptotic theory
fossil	Palaeoecological and Palaeogeographical Analysis Tools
fourPNO	Bayesian 4 Parameter Item Response Model
FourScores	FourScores - A game for two players
fpc	Flexible Procedures for Clustering
fpca	Restricted MLE for Functional Principal Components Analysis
fpCompare	Reliable Comparison of Floating Point Numbers
FPDclustering	PD-Clustering and Factor PD-Clustering
fPortfolio	Rmetrics - Portfolio Selection and Optimization
fpow	Computing the noncentrality parameter of the noncentral F distribution
fpp	Data for "Forecasting: principles and practice"
fptdApprox	Approximation of First-Passage-Time Densities for Diffusion Processes
fracdiff	Fractionally differenced ARIMA aka ARFIMA(p,d,q) models
fracprolif	Fraction Proliferation via a Quiescent Growth Model
fractal	Fractal Time Series Modeling and Analysis
fractaldim	Estimation of fractal dimensions
FractalParameterEstimation	Estimation of Parameters p and q for Randomized Sierpinski Carpet for [p-p-p-q]-Model
fractalrock	Generate fractal time series with non-normal returns distribution
FRACTION	Numeric number into fraction
fractional	Vulgar Fractions in R
Fragman	Fragment Analysis in R
frailtyHL	Frailty Models via H-likelihood
frailtypack	General Frailty Models: Shared, Joint and Nested Frailty Models with Prediction
frailtySurv	General Semiparametric Shared Frailty Model
frair	Functional response analysis in R
Frames2	Estimation in Dual Frame Surveys
franc	Detect the Language of Text
FRAPO	Financial Risk Modelling and Portfolio Optimisation with R
FRB	Fast and Robust Bootstrap
frbs	Fuzzy Rule-Based Systems for Classification and Regression Tasks
FRCC	Fast Regularized Canonical Correlation Analysis
freeknotsplines	Free-Knot Splines
FreeSortR	Free Sorting data analysis
freestats	Statistical algorithms used in common data mining course
FREGAT	Family REGional Association Tests
fRegression	Rmetrics - Regression Based Decision and Prediction
FREQ	FREQ: Estimate population size from capture frequencies
freqdom	Frequency Domain Analysis for Multivariate Time Series
freqMAP	Frequency Moving Average Plots (MAP) of Multinomial Data by a Continuous Covariate
freqparcoord	Novel Methods for Parallel Coordinates
FreqProf	Frequency Profiles Computing and Plotting
freqweights	Working with Frequency Tables
FRESA.CAD	Feature Selection Algorithms for Computer Aided Diagnosis
FrF2	Fractional Factorial designs with 2-level factors
FrF2.catlg128	Catalogues of resolution IV 128 run 2-level fractional factorials up to 33 factors that do have 5-letter words
frm	Regression Analysis of Fractional Responses
frmhet	Regression Analysis of Fractional Responses Under Unobserved Heterogeneity
frmpd	Regression Analysis of Panel Fractional Responses
frmqa	The Generalized Hyperbolic Distribution, Related Distributions and Their Applications in Finance
frontier	Stochastic Frontier Analysis
frontiles	Partial Frontier Efficiency Analysis
frt	Full Randomization Test
FSA	Functions for Simple Fisheries Stock Assessment Methods
FSAdata	Data to Support Fish Stock Assessment (FSA) Package
fscaret	Automated Feature Selection from 'caret'
FSelector	Selecting attributes
fsia	Import and Analysis of OMR Data from FormScanner
FSInteract	Fast Searches for Interactions
fslr	Wrapper Functions for FSL (FMRIB Software Library) from Functional MRI of the Brain (FMRIB)
fso	Fuzzy Set Ordination
fSRM	Social Relations Analyses with Roles ("Family SRM")
FTICRMS	Programs for Analyzing Fourier Transform-Ion Cyclotron Resonance Mass Spectrometry Data
ftnonpar	Features and Strings for Nonparametric Regression
fTrading	Technical Trading Analysis
fts	R interface to tslib (a time series library in c++)
ftsa	Functional Time Series Analysis
ftsspec	Spectral Density Estimation and Comparison for Functional Time Series
fueleconomy	EPA fuel economy data
fugeR	FUzzy GEnetic, a machine learning algorithm to construct prediction model based on fuzzy logic
fullfact	Full Factorial Breeding Analysis
fulltext	Full Text of 'Scholarly' Articles Across Many Data Sources
fun	Use R for Fun
FunChisq	Chi-Square and Exact Tests for Non-Parametric Functional Dependencies
FunCluster	Functional Profiling of Microarray Expression Data
Funclustering	A package for functional data clustering
FuncMap	Hive Plots of R Package Function Calls
functional	Curry, Compose, and other higher-order functions
FunctionalNetworks	An algorithm for gene and gene set network inference
functools	Functional Programming in R
funcy	Functional Clustering Algorithms
funFEM	Clustering in the Discriminative Functional Subspace
fungible	Fungible Coefficients and Monte Carlo Functions
funHDDC	Model-based clustering in group-specific functional subspaces
fUnitRoots	Trends and Unit Roots
funModeling	Learning Data Cleaning, Visual Analysis and Model Performance
funr	Simple Utility Providing Terminal Access to all R Functions
funreg	funreg (Functional Regression for Irregularly Timed Data)
funtimes	Functions for Time Series Analysis
FusedPCA	Community Detection via Fused Principal Component Analysis
futile.any	A Tiny Utility Providing Polymorphic Operations
futile.logger	A Logging Utility for R
futile.matrix	Random matrix generation and manipulation
futile.options	Futile options management
futile.paradigm	A framework for working in a functional programming paradigm in R
future	A Future API for R
futureheatwaves	Find, Characterize, and Explore Heat Waves in Climate Projections
fuzzyFDR	Exact calculation of fuzzy decision rules for multiple testing
fuzzyforest	Fuzzy Forests
FuzzyLP	Fuzzy Linear Programming
FuzzyNumbers	Tools to Deal with Fuzzy Numbers
fuzzyRankTests	Fuzzy Rank Tests and Confidence Intervals
FuzzyStatProb	Fuzzy Stationary Probabilities from a Sequence of Observations of an Unknown Markov Chain
FuzzyToolkitUoN	Type 1 Fuzzy Logic Toolkit
fwdmsa	Forward search for Mokken scale analysis
FWDselect	Selecting Variables in Regression Models
fwi.fbp	Fire Weather Index System and Fire Behaviour Prediction System Calculations
fwsim	Fisher-Wright Population Simulation
fxregime	Exchange Rate Regime Analysis
0	
G1DBN	A package performing Dynamic Bayesian Network inference
G2Sd	Grain-Size Statistics and Description of Sediment
GA	Genetic Algorithms
GA4Stratification	A genetic algorithm approach to determine stratum boundaries and sample sizes of each stratum in stratified sampling
GAabbreviate	Abbreviating Items Measures using Genetic Algorithms
GABi	Framework for Generalized Subspace Pattern Mining
GAD	GAD: Analysis of variance from general principles
GaDiFPT	First Passage Time Simulation for Gaussian Diffusion Processes
gains	Gains Table Package
GAIPE	Graphical Extension with Accuracy in Parameter Estimation (GAIPE)
galts	Genetic algorithms and C-steps based LTS (Least Trimmed Squares) estimation
gam	Generalized Additive Models
gamair	Data for "GAMs: An Introduction with R"
gambin	Fit the GamBin Model to Species Abundance Distributions
GAMBoost	Generalized linear and additive models by likelihood based boosting
gamboostLSS	Boosting Methods for 'GAMLSS'
gamboostMSM	Estimating multistate models using gamboost()
gamclass	Functions and Data for a Course on Modern Regression and Classification
GAMens	Applies GAMbag, GAMrsm and GAMens Ensemble Classifiers for Binary Classification
games	Statistical Estimation of Game-Theoretic Models
GameTheory	Cooperative Game Theory
gamlr	Gamma Lasso Regression
gamlss	Generalised Additive Models for Location Scale and Shape
gamlss.add	Extra Additive Terms for GAMLSS Models
gamlss.cens	Fitting an Interval Response Variable Using gamlss.family Distributions
gamlss.data	GAMLSS Data
gamlss.demo	Demos for GAMLSS
gamlss.dist	Distributions to be Used for GAMLSS Modelling
gamlss.mx	Fitting Mixture Distributions with GAMLSS
gamlss.nl	Fitting non linear parametric GAMLSS models
gamlss.spatial	Spatial Terms in GAMLSS Models
gamlss.tr	Generating and fitting truncated (gamlss.family) distributions
gamlss.util	GAMLSS Utilities
gamm4	Generalized additive mixed models using mgcv and lme4
Gammareg	classic gamma regression: joint modeling of mean and shape parameters
gammSlice	Generalized additive mixed model analysis via slice sampling
gamsel	Fit Regularization Path for Generalized Additive Models
GANPA	Gene Association Network-based Pathway Analysis
GANPAdata	The GANPA Datasets Package
gaoptim	Genetic Algorithm optimization for real-based and permutation-based problems
gap	Genetic Analysis Package
gapmap	Functions for Drawing Gapped Cluster Heatmap with ggplot2
gapminder	Data from Gapminder
GAR	Authorize and Request Google Analytics Data
gaselect	Genetic Algorithm (GA) for Variable Selection from High-Dimensional Data
gaston	Genetic Data Manipulation (Quality Control, GRM and LD Computations, PCA), Linear Mixed Models (AIREML Algorithm), Association Testing
gaussDiff	Difference measures for multivariate Gaussian probability density functions
gaussquad	Collection of functions for Gaussian quadrature
gazepath	Gazepath Transforms Eye-Tracking Data into Fixations and Saccades
gb	Generalize Lambda Distribution and Generalized Bootstrapping
GB2	Generalized Beta Distribution of the Second Kind: Properties, Likelihood, Estimation
gbm	Generalized Boosted Regression Models
gbm2sas	Convert GBM Object Trees to SAS Code
gbRd	Utilities for processing Rd objects and files
gbutils	Simulation of Real and Complex Numbers and Small Programming Utilities
GCAI.bias	Guided Correction Approach for Inherited bias (GCAI.bias)
gCat	Graph-based two-sample tests for categorical data
gcbd	GPU/CPU Benchmarking in Debian-based systems
GCD	Global Charcoal Database
gcdnet	LASSO and (adaptive) Elastic-Net penalized least squares, logistic regression, HHSVM and squared hinge loss SVM using a fast GCD algorithm
gcerisk	Generalized Competing Event Model
gclus	Clustering Graphics
gcmr	Gaussian Copula Marginal Regression
gconcord	Concord method for Graphical Model Selection
gcookbook	Data for "R Graphics Cookbook"
GCPM	Generalized Credit Portfolio Model
GDAdata	Datasets for the Book Graphical Data Analysis with R
gdalUtils	Wrappers for the Geospatial Data Abstraction Library (GDAL) Utilities
gdata	Various R Programming Tools for Data Manipulation
g.data	Delayed-Data Packages
GDAtools	A toolbox for the analysis of categorical data in social sciences, and especially Geometric Data Analysis
GDELTtools	Download, slice, and normalize GDELT data
gdimap	Generalized Diffusion Magnetic Resonance Imaging
gdistance	Distances and Routes on Geographical Grids
gdm	Functions for Generalized Dissimilarity Modeling
gdtools	Utilities for Graphical Rendering
gear	Geostatistical Analysis in R
gee	Generalized Estimation Equation Solver
geeM	Solve Generalized Estimating Equations
geepack	Generalized Estimating Equation Package
geesmv	Modified Variance Estimators for Generalized Estimating Equations
geigen	Calculate Generalized Eigenvalues, the Generalized Schur Decomposition and the Generalized Singular Value Decomposition of a Matrix Pair with Lapack
geiger	Analysis of Evolutionary Diversification
gelnet	Generalized Elastic Nets
gems	Generalized Multistate Simulation Model
gemtc	Network Meta-Analysis Using Bayesian Methods
gemtc.jar	GeMTC Java binary
GenABEL	genome-wide SNP association analysis
GenABEL.data	Package contains data which is used by GenABEL example and test functions
genalg	R Based Genetic Algorithm
genasis	Global ENvironmental ASsessment Information System (GENASIS) computational tools
GenBinomApps	Clopper-Pearson Confidence Interval and Generalized Binomial Distribution
GenCAT	Genetic Class Association Testing
gendata	Generate and Modify Synthetic Datasets
gender	Predict Gender from Names Using Historical Data
genderizeR	Gender Prediction Based on First Names
gendist	Generated Probability Distribution Models
GENEAread	Package For Reading Binary files
GeneCycle	Identification of Periodically Expressed Genes
GeneF	Package for Generalized F-statistics
GeneFeST	Bayesian calculation of gene-specific FST from genomic SNP data
Geneland	Detection of structure from multilocus genetic data
geneListPie	Profiling a gene list into GOslim or KEGG function pie
GeneNet	Modeling and Inferring Gene Networks
geneNetBP	Belief Propagation in Genotype-Phenotype Networks
genepi	Genetic Epidemiology Design and Inference
GeneralizedHyperbolic	The Generalized Hyperbolic Distribution
GeneralOaxaca	Blinder-Oaxaca Decomposition for Generalized Linear Model
generator	Generate Data Containing Fake Personally Identifiable Information
GeneReg	Construct time delay gene regulatory network
geneSignatureFinder	A Gene-signatures finder tools
geneSLOPE	Genome-Wide Association Study with SLOPE
genetics	Population Genetics
GeneticSubsetter	Identify Favorable Subsets of Germplasm Collections
GeneticTools	Collection of Genetic Data Analysis Tools
GeNetIt	Spatial Graph-Theoretic Genetic Gravity Modelling
GenForImp	The Forward Imputation: A Sequential Distance-Based Approach for Imputing Missing Data
genie	A New, Fast, and Outlier Resistant Hierarchical Clustering Algorithm
GenKern	Functions for generating and manipulating binned kernel density estimates
genlasso	Path algorithm for generalized lasso problems
GENLIB	Genealogical Data Analysis
genMOSS	Functions for the Bayesian Analysis of GWAS Data
genMOSSplus	Application of MOSS algorithm to genome-wide association study (GWAS)
genoPlotR	Plot Publication-Grade Gene and Genome Maps
GenOrd	Simulation of Discrete Random Variables with Given Correlation Matrix and Marginal Distributions
genpathmox	Generalized PATHMOX Algorithm for PLS-PM, LS and LAD Regression
genridge	Generalized Ridge Trace Plots for Ridge Regression
GenSA	R Functions for Generalized Simulated Annealing
gensemble	generalized ensemble methods
genSurv	Generating Multi-State Survival Data
GenWin	Spline Based Window Boundaries for Genomic Analyses
geo	Draw and Annotate Maps, Especially Charts of the North Atlantic
geoaxe	Split 'Geospatial' Objects into Pieces
geoBayes	Analysis of Geostatistical Data using Bayes and Empirical Bayes Methods
GeoBoxplot	Geographic Box Plot
geocodeHERE	Wrapper for Nokia's HERE Geocoding API
geoCount	Analysis and Modeling for Geostatistical Count Data
GeoDE	A geometrical Approach to Differential expression and gene-set enrichment
geoelectrics	3D-Visualization of Geoelectric Resistivity Measurement Profiles
geofd	Spatial Prediction for Function Value Data
GeoGenetix	Quantification of the effect of geographic versus environmental isolation on genetic differentiation
geojsonio	Convert Data from and to 'geoJSON' or 'topoJSON'
geoknife	Web-Processing of Large Gridded Datasets
GeoLight	Analysis of Light Based Geolocator Data
GEOmap	Topographic and Geologic Mapping
geomapdata	Data for topographic and Geologic Mapping
geometry	Mesh Generation and Surface Tesselation
geomnet	Network Visualization in the 'ggplot2' Framework
geomorph	Geometric Morphometric Analyses of 2D/3D Landmark Data
geonames	Interface to www.geonames.org web service
geophys	Geophysics, Continuum Mechanics, Mogi Models, Gravity
geoR	Analysis of Geostatistical Data
geoRglm	A Package for Generalised Linear Spatial Models
georob	Robust Geostatistical Analysis of Spatial Data
geoscale	Geological Time Scale Plotting
geospacom	Facilitate Generating of Distance Matrices Used in Package 'spacom' and Plotting Data on Maps
geosphere	Spherical Trigonometry
geospt	Geostatistical Analysis and Design of Optimal Spatial Sampling Networks
geosptdb	Spatio-Temporal; Inverse Distance Weighting and Radial Basis Functions with Distance-Based Regression
geostatsp	Geostatistical Modelling with Likelihood and Bayes
geotech	Geotechnical Engineering
geotools	Geo tools
geotopbricks	An R Plug-in for the Distributed Hydrological Model GEOtop
GeoXp	Interactive exploratory spatial data analysis
geozoo	Zoo of Geometric Objects
GERGM	Estimation and Fit Diagnostics for Generalized Exponential Random Graph Models
gesca	Generalized Structured Component Analysis (GSCA)
gesis	R Client for GESIS Data Catalogue (DBK)
GESTr	Gene Expression State Transformation
getMet	Get Meteorological Data for Hydrologic Models
getopt	C-like getopt behavior
GetoptLong	Parsing Command-Line Arguments and Variable Interpolation
getPass	Masked User Input
GetR	GetR: Calculate Guttman error trees in R
gets	General-to-Specific (GETS) Modelling and Indicator Saturation Methods
GetTDData	Get Data for Brazilian Bonds (Tesouro Direto)
gettingtothebottom	Learning Optimization and Machine Learning for Statistics
GEVcdn	GEV Conditional Density Estimation Network
GEVStableGarch	ARMA-GARCH/APARCH Models with GEV and Stable Distributions
GExMap	A visual, intuitive, easy to use software giving access to a new type of information buried into your microarray data
gfcanalysis	Tools for Working with Hansen et al. Global Forest Change Dataset
GFD	Tests for General Factorial Designs
GGally	Extension to ggplot2
ggalt	Extra Coordinate Systems, Geoms and Statistical Transformations for 'ggplot2'
ggbeeswarm	Categorical Scatter (Violin Point) Plots
ggcorrplot	Visualization of a Correlation Matrix using 'ggplot2'
ggdendro	Create Dendrograms and Tree Diagrams Using 'ggplot2'
gge	Genotype Plus Genotype-by-Environment Biplots
GGEBiplotGUI	Interactive GGE Biplots in R
ggenealogy	Visualization Tools for Genealogical Data
ggExtra	Add Marginal Histograms to 'ggplot2', and More 'ggplot2' Enhancements
ggfortify	Data Visualization Tools for Statistical Analysis Results
GGIR	Raw Accelerometer Data Analysis
ggiraph	Make 'ggplot2' Graphics Interactive Using 'htmlwidgets'
gglasso	Group Lasso Penalized Learning Using A Unified BMD Algorithm
ggm	Functions for graphical Markov models
ggmap	Spatial Visualization with ggplot2
ggmcmc	Tools for Analyzing MCMC Simulations from Bayesian Inference
GGMridge	Gaussian Graphical Models Using Ridge Penalty Followed by Thresholding and Reestimation
GGMselect	Gaussian Graphs Models Selection
ggnetwork	Geometries to Plot Networks with 'ggplot2'
ggparallel	Variations of Parallel Coordinate Plots for Categorical Data
ggplot2	An Implementation of the Grammar of Graphics
ggplot2movies	Movies Data
ggpmisc	Miscellaneous Extensions to 'ggplot2'
ggRandomForests	Visually Exploring Random Forests
ggraptR	Allows Interactive Visualization of Data Through a Web Browser GUI
ggrepel	Repulsive Text and Label Geoms for 'ggplot2'
ggROC	package for roc curve plot with ggplot2
ggseas	'stats' for Seasonal Adjustment on the Fly with 'ggplot2'
ggsn	North Symbols and Scale Bars for Maps Created with 'ggplot2' or 'ggmap'
ggspectra	Extensions to 'ggplot2' for Radiation Spectra
ggswissmaps	Offers Various Swiss Maps as ggplot2 Objects
ggtern	An Extension to 'ggplot2', for the Creation of Ternary Diagrams
ggThemeAssist	Add-in to Customize 'ggplot2' Themes
ggthemes	Extra Themes, Scales and Geoms for 'ggplot2'
ggvis	Interactive Grammar of Graphics
ghit	Lightweight GitHub Package Installer
GHQp	Gauss Hermite Quadrature with pruning
ghyp	A package on the generalized hyperbolic distribution and its special cases
GiANT	Gene Set Uncertainty in Enrichment Analysis
GibbsACOV	Gibbs Sampler for One-Way Mixed-Effects ANOVA and ANCOVA Models
gibbs.met	Naive Gibbs Sampling with Metropolis Steps
GIGrvg	Random Variate Generator for the GIG Distribution
GillespieSSA	Gillespie's Stochastic Simulation Algorithm (SSA)
gimme	Group Iterative Multiple Model Estimation
gimms	Download and Process GIMMS NDVI3g Data
GiniWegNeg	Computing the Gini Coefficient for Weighted and Negative Attributes
gIPFrm	Generalized Iterative Proportional Fitting for Relational Models
GiRaF	Gibbs Random Fields Analysis
giRaph	The giRaph package for graph representation in R
GISTools	Some further GIS capabilities for R
gistr	Work with 'GitHub' 'Gists'
git2r	Provides Access to Git Repositories
gitlabr	Access to the Gitlab API
gitter	Quantification of Pinned Microbial Cultures
Giza	Constructing panels of population pyramid plots based on lattice
gjam	Generalized Joint Attribute Modeling
gkmSVM	Gapped-Kmer Support Vector Machine
glamlasso	Lasso Penalization in Large Scale Generalized Linear Array Models
glarma	Generalized Linear Autoregressive Moving Average Models
glasso	Graphical lasso- estimation of Gaussian graphical models
glba	General Linear Ballistic Accumulator Models
glcm	Calculate Textures from Grey-Level Co-Occurrence Matrices (GLCMs)
gld	Estimation and Use of the Generalised (Tukey) Lambda Distribution
GLDEX	Fitting Single and Mixture of Generalised Lambda Distributions (RS and FMKL) using Various Methods
gldist	An Asymmetry-Steepness Parameterization of the Generalized Lambda Distribution
GLDreg	Fit GLD Regression Model and GLD Quantile Regression Model to Empirical Data
glinternet	Learning Interactions via Hierarchical Group-Lasso Regularization
gllm	Generalised log-linear model
glm2	Fitting Generalized Linear Models
glmc	Fitting Generalized Linear Models Subject to Constraints
glm.ddR	Distributed 'glm' for Big Data using 'ddR' API
glmdm	R Code for Simulation of GLMDM
glmgraph	Graph-Constrained Regularization for Sparse Generalized Linear Models
glmlep	Fit GLM with LEP-based penalized maximum likelihood
glmm	Generalized Linear Mixed Models via Monte Carlo Likelihood Approximation
glmmBUGS	Generalised Linear Mixed Models and Spatial Models with WinBUGS, BRugs, or OpenBUGS
glmmGS	Gauss-Seidel Generalized Linear Mixed Model solver
glmmLasso	Variable Selection for Generalized Linear Mixed Models by L1-Penalized Estimation
glmmML	Generalized linear models with clustering
GLMMRR	Generalized Linear Mixed Model (GLMM) for Binary Randomized Response Data
glmmsr	Fit a Generalized Linear Mixed Model
glmnet	Lasso and Elastic-Net Regularized Generalized Linear Models
glmnetcr	Fit a penalized constrained continuation ratio model for predicting an ordinal response
glmpath	L1 Regularization Path for Generalized Linear Models and Cox Proportional Hazards Model
glmpathcr	Fit a penalized continuation ratio model for predicting an ordinal response
glmulti	Model selection and multimodel inference made easy
glmvsd	Variable Selection Deviation Measures and Instability Tests for High-Dimensional Generalized Linear Models
glmx	Generalized Linear Models Extended
globalboosttest	Testing the additional predictive value of high-dimensional data
GlobalDeviance	Global Deviance Permutation Tests
GlobalFit	Bi-Level Optimization of Metabolic Network Models
globalGSA	Global Gene-Set Analysis for Association Studies
GlobalOptions	Generate Functions to Get or Set Global Options
globalOptTests	Objective functions for benchmarking the performance of global optimization algorithms
globals	Identify Global Objects in R Expressions
globe	Plot 2D and 3D Views of the Earth, Including Major Coastline
glogis	Fitting and Testing Generalized Logistic Distributions
glpkAPI	R Interface to C API of GLPK
glrt	Generalized Logrank Tests for Interval-censored Failure Time Data
GLSME	Generalized Least Squares with Measurement Error
glycanr	Tools for Analysing N-Glycan Data
gmailr	Access the Gmail RESTful API
gmapsdistance	Distance and Travel Time Between Two Points from Google Maps
gmatrix	GPU Computing in R
GMCM	Fast Estimation of Gaussian Mixture Copula Models
gMCP	Graph Based Multiple Comparison Procedures
GMD	Generalized Minimum Distance of distributions
GMDH	Short Term Forecasting via GMDH-Type Neural Network Algorithms
Gmedian	Geometric Median, k-Median Clustering and Robust Median PCA
gmeta	Meta-Analysis via a Unified Framework of Confidence Distribution
Gmisc	Descriptive Statistics, Transition Plots, and More
gmm	Generalized Method of Moments and Generalized Empirical Likelihood
GMMBoost	Likelihood-based Boosting for Generalized mixed models
gmnl	Multinomial Logit Models with Random Parameters
gmodels	Various R Programming Tools for Model Fitting
gmp	Multiple Precision Arithmetic
gmt	Interface between GMT Map-Making Software and R
gmum.r	GMUM Machine Learning Group Package
gmwm	Generalized Method of Wavelet Moments
gMWT	Generalized Mann-Whitney Type Tests
GNE	Computation of Generalized Nash Equilibria
gnm	Generalized Nonlinear Models
gnmf	Generalized Non-negative Matrix Factorization Based on Renyi Divergence
gnumeric	Read Data from Files Readable by Gnumeric
goalprog	Weighted and lexicographical goal programming and optimization
gof	Model-diagnostics based on cumulative residuals
gofCopula	Goodness-of-Fit Tests for Copulae
GoFKernel	Testing Goodness-of-Fit with the Kernel Density Estimator
goft	Tests of Fit for some Probability Distributions
goftest	Classical Goodness-of-Fit Tests for Univariate Distributions
GOGANPA	GO-Functional-Network-based Gene-Set-Analysis
gogarch	Generalized Orthogonal GARCH (GO-GARCH) models
GoodmanKruskal	Association Analysis for Categorical Variables
googleAuthR	Easy Authentication with Google OAuth2 APIs
googleformr	Collect Data Programmatically by POST Methods to Google Forms
googlePublicData	Working with Google Public Data Explorer DSPL Metadata Files
googlesheets	Manage Google Spreadsheets from R
googleVis	R Interface to Google Charts
GOplot	Visualization of Functional Analysis Data
GORCure	Fit Generalized Odds Rate Mixture Cure Model with Interval Censored Data
goric	Generalized Order-Restricted Information Criterion
Goslate	Goslate Interface
govStatJPN	functions to get public survey data in Japan
gpairs	gpairs: The Generalized Pairs Plot
GPareto	Gaussian Processes for Pareto Front Estimation and Optimization
GPArotation	GPA Factor Rotation
GPC	Generalized Polynomial Chaos
gPCA	Batch Effect Detection via Guided Principal Components Analysis
gpclib	General Polygon Clipping Library for R
GPCSIV	GPCSIV, Generalized Principal Component of Symbolic Interval variables
gpDDE	General Profiling Method for Delay Differential Equation
gPdtest	Bootstrap goodness-of-fit test for the generalized Pareto distribution
GPFDA	Apply Gaussian Process in Functional data analysis
GPfit	Gaussian Processes Modeling
gpk	100 Data Sets for Statistics Education
gplm	Generalized partial linear models (GPLM)
gplots	Various R Programming Tools for Plotting Data
GPLTR	Generalized Partially Linear Tree-Based Regression Model
gpmap	Analysing and plotting genotype-phenotype maps
gpr	A Minimalistic package to apply Gaussian Process in R
gProfileR	Interface to the 'g:Profiler' Toolkit
GPseq	gpseq: Using the generalized Poisson distribution to model sequence read counts from high throughput sequencing experiments
gptk	Gaussian Processes Tool-Kit
gpuR	GPU Functions for R Objects
gputools	A Few GPU Enabled Functions
GPvam	Maximum Likelihood Estimation of Multiple Membership Mixed Models Used in Value-Added Modeling
gquad	G Quadruplex Motif Prediction Tool
Grace	Graph-Constrained Estimation and Hypothesis Tests
grade	Binary Grading functions for R
GRaF	Species distribution modelling using latent Gaussian random fields
gRain	Graphical Independence Networks
gramEvol	Grammatical Evolution for R
GrammR	Graphical Representation and Modeling of Metagenomic Reads
granova	Graphical Analysis of Variance
granovaGG	Graphical Analysis of Variance Using ggplot2
gRapfa	Acyclic Probabilistic Finite Automata
gRapHD	Efficient selection of undirected graphical models for high-dimensional datasets
GrapheR	A Multi-Platform GUI for Drawing Customizable Graphs in R
graphicalVAR	Graphical VAR for Experience Sampling Data
graphicsQC	Quality Control for Graphics in R
graphscan	Cluster Detection with Hypothesis Free Scan Statistic
graphTweets	Visualise Twitter Interactions
GrassmannOptim	Grassmann Manifold Optimization
graticule	Meridional and Parallel Lines for Maps
gRbase	A Package for Graphical Modelling in R
gRc	Inference in Graphical Gaussian Models with Edge and Vertex Symmetries
Greg	Regression Helper Functions
greport	Graphical Reporting for Clinical Trials
greyzoneSurv	Fit a Grey-Zone Model with Survival Data
Grid2Polygons	Convert Spatial Grids to Polygons
gridBase	Integration of base and grid graphics
gridDebug	Debugging 'grid' Graphics
gridExtra	Miscellaneous Functions for "Grid" Graphics
gridGraphics	Redraw Base Graphics Using 'grid' Graphics
gridGraphviz	Drawing Graphs with 'grid'
gridSVG	Export 'grid' Graphics as SVG
GriegSmith	Uses Grieg-Smith method on 2 dimentional spatial data
gRim	Graphical Interaction Models
grImport	Importing Vector Graphics
grnn	General regression neural network
groc	Generalized Regression on Orthogonal Components
grofit	The package was developed to fit fit many growth curves obtained under different conditions
gromovlab	Gromov-Hausdorff Type Distances for Labeled Metric Spaces
grouped	Regression Analysis of Grouped and Coarse Data
groupRemMap	Regularized Multivariate Regression for Identifying Master Predictors Using the GroupRemMap Penalty
GroupSeq	A GUI-Based Program to Compute Probabilities Regarding Group Sequential Designs
GroupTest	Multiple Testing Procedure for Grouped Hypotheses
growcurves	Bayesian Semi and Nonparametric Growth Curve Models that Additionally Include Multiple Membership Random Effects
growfunctions	Bayesian Non-Parametric Dependent Models for Time-Indexed Functional Data
growthcurver	Simple Metrics to Summarize Growth Curves
growthmodels	Nonlinear Growth Models
growthrate	Bayesian reconstruction of growth velocity
growthrates	Estimate Growth Rates from Experimental Data
grplasso	Fitting user specified models with Group Lasso penalty
grppenalty	Concave 1-norm and 2-norm group penalty in linear and logistic regression
grpreg	Regularization Paths for Regression Models with Grouped Covariates
grpregOverlap	Penalized Regression Models with Overlapping Grouped Covariates
grpss	Group Screening and Selection
grr	Alternative Implementations of Base R Functions
grt	General Recognition Theory
GRTo	Tools for the Analysis of Gutenberg-Richter Distributions of Earthquake Magnitudes
GSA	Gene set analysis
GSAgm	Gene Set Analysis using the Gamma Method
gsalib	Utility Functions For GATK
gsarima	Two functions for Generalized SARIMA time series simulation
gsbDesign	Group Sequential Bayes Design
gsDesign	Group Sequential Design
GSE	Robust Estimation in the Presence of Cellwise and Casewise Contamination and Missing Data
gsEasy	Gene Set Enrichment Analysis in R
gSeg	Graph-Based Change-Point Detection (g-Segmentation)
gSEM	Semi-Supervised Generalized Structural Equation Modeling
gset	Group Sequential Design in Equivalence Studies
gsg	Calculation of selection coefficients
gsheet	Download Google Sheets Using Just the URL
GSIF	Global Soil Information Facilities
gskat	GEE_KM
gsl	wrapper for the Gnu Scientific Library
GSM	Gamma Shape Mixture
gsmoothr	Smoothing tools
gss	General Smoothing Splines
gsscopu	Copula Density and 2-D Hazard Estimation using Smoothing Splines
GSSE	Genotype-Specific Survival Estimation
gstat	Spatial and Spatio-Temporal Geostatistical Modelling, Prediction and Simulation
gsubfn	Utilities for strings and function arguments
gsw	Gibbs Sea Water Functions
GsymPoint	Estimation of the Generalized Symmetry Point, an Optimal Cutpoint in Continuous Diagnostic Tests
gtable	Arrange 'Grobs' in Tables
gtcorr	Calculate efficiencies of group testing algorithms with correlated responses
gte	Generalized Turnbull's Estimator
gtools	Various R Programming Tools
gtop	Game-Theoretically OPtimal (GTOP) Reconciliation Method
gtrendsR	R Functions to Perform and Display Google Trends Queries
gtx	Genetics ToolboX
GuardianR	The Guardian API Wrapper
Guerry	Maps, data and methods related to Guerry (1833) "Moral Statistics of France"
guess	Adjust Estimates of Learning for Guessing
GUIDE	GUI for DErivatives in R
GUILDS	Implementation of Sampling Formulas for the Unified Neutral Model of Biodiversity and Biogeography, with or without Guild Structure
GUIProfiler	Graphical User Interface for Rprof()
gumbel	The Gumbel-Hougaard Copula
GUniFrac	Generalized UniFrac distances
gunsales	Statistical Analysis of Monthly Background Checks of Gun Purchases
GUTS	Fast Calculation of the Likelihood of a Stochastic Survival Model
gvc	Global Value Chains Tools
gvcm.cat	Regularized Categorical Effects/Categorical Effect Modifiers/Continuous/Smooth Effects in GLMs
gvlma	Global Validation of Linear Models Assumptions
GWAF	Genome-Wide Association/Interaction Analysis and Rare Variant Analysis with Family Data
GWASExactHW	Exact Hardy-Weinburg testing for Genome Wide Association Studies
gwerAM	Controlling the genome-wide type I error rate in association mapping experiments
GWG	Calculation of probabilities for inadequate and excessive gestational weight gain
gWidgets	gWidgets API for building toolkit-independent, interactive GUIs
gWidgets2	Rewrite of gWidgets API for Simplified GUI Construction
gWidgets2RGtk2	Implementation of gWidgets2 for the RGtk2 Package
gWidgets2tcltk	Toolkit Implementation of gWidgets2 for tcltk
gWidgetsRGtk2	Toolkit implementation of gWidgets for RGtk2
gWidgetstcltk	Toolkit implementation of gWidgets for tcltk package
GWLelast	Geographically Weighted Logistic Elastic Net Regression
GWmodel	Geographically-Weighted Models
GWRM	Generalized Waring Regression Model for Count Data
gwrr	Fits geographically weighted regression models with diagnostic tools
GWsignif	Genome-wide significance for whole genome sequencing studies
GxM	Maximum Likelihood Estimation for Gene-by-Measured Environment Interaction Models
gyriq	Kinship-Adjusted Survival SNP-Set Analysis
0	
h2o	R Interface for H2O
h5	Interface to the 'HDF5' Library
haarfisz	Software to perform Haar Fisz transforms
HAC	Estimation, Simulation and Visualization of Hierarchical Archimedean Copulae (HAC)
HadoopStreaming	Utilities for using R scripts in Hadoop streaming
hail	Read HYDRA Rainfall Data
hamlet	Hierarchical Optimal Matching and Machine Learning Toolbox
HandTill2001	Multiple Class Area under ROC Curve
Hankel	Univariate non-parametric two-sample test based on empirical Hankel transforms
hapassoc	Inference of Trait Associations with SNP Haplotypes and Other Attributes using the EM Algorithm
HapEstXXR	Multi-Locus Stepwise Regression
HAPim	HapIM
Haplin	Analyzing Case-Parent Triad and/or Case-Control Data with SNP Haplotypes
haplo.ccs	Estimate Haplotype Relative Risks in Case-Control Data
HaploSim	Functions to simulate haplotypes
haplo.stats	Statistical Analysis of Haplotypes with Traits and Covariates when Linkage Phase is Ambiguous
haplotypes	Haplotype Inference and Statistical Analysis of Genetic Variation
HAP.ROR	Recursive Organizer (ROR)
HardyWeinberg	Graphical Tests for Hardy-Weinberg Equilibrium
HarmonicRegression	Harmonic Regression to One or more Time Series
harvestr	A Parallel Simulation Framework
Harvest.Tree	Harvest the Classification Tree
hash	Full feature implementation of hash/associated arrays/dictionaries
hashFunction	A collection of non-cryptographic hash functions
hashids	Generate Short Unique YouTube-Like IDs (Hashes) from Integers
hashr	Hash R Objects to Integers Fast
hasseDiagram	Drawing Hasse diagram
haven	Import SPSS, Stata and SAS Files
hawkes	Hawkes process simulation and calibration toolkit
hazus	Damage functions from FEMA's HAZUS software for use in modeling financial losses from natural disasters
hBayesDM	Hierarchical Bayesian Modeling of Decision-Making Tasks
HBglm	Hierarchical Bayesian Regression for GLMs
hbim	Hill/Bliss Independence Model for Combination Vaccines
hbm	Hierarchical Block Matrix Analysis
hbmem	Hierarchical Bayesian Analysis of Recognition Memory
hbsae	Hierarchical Bayesian Small Area Estimation
HBSTM	Hierarchical Bayesian Space-Time models for Gaussian space-time data
hcc	Hidden correlation check
hcci	Interval estimation for the parameters of linear models with heteroskedasticity (Wild Bootstrap)
hcp	Change Point Estimation for Regression with Heteroscedastic Data
hda	Heteroscedastic Discriminant Analysis
HDclassif	High Dimensional Supervised Classification and Clustering
hddplot	Use known groups in high-dimensional data to derive scores for plots
hddtools	Hydrological Data Discovery Tools
hdeco	Hierarchical DECOmposition of Entropy for Categorical Map Comparisons
HDGLM	Tests for High Dimensional Generalized Linear Models
hdi	High-Dimensional Inference
HDInterval	Highest (Posterior) Density Intervals
hdlm	Fitting High Dimensional Linear Models
hdm	High-Dimensional Metrics
HDMD	Statistical Analysis Tools for High Dimension Molecular Data (HDMD)
hdnom	Nomograms for High-Dimensional Cox Models
HDPenReg	High-Dimensional Penalized Regression
hdr	Interface to the UNDR Human Development Report API
hdrcde	Highest density regions and conditional density estimation
HDtweedie	The Lasso for the Tweedie's Compound Poisson Model Using an IRLS-BMD Algorithm
HEAT	Health Effects of Air Pollution and Temperature (HEAT)
heatex	Heat exchange calculations during physical activity
heatmap3	An Improved Heatmap Package
heatmapFit	Fit statistic for binary dependent variable models
heatmap.plus	Heatmap with more sensible behavior
heavy	Robust Estimation Using Heavy-Tailed Distributions
heemod	Models for Health Economic Evaluation
hellno	Providing 'stringsAsFactors=FALSE' Variants of 'data.frame()' and 'as.data.frame()'
helloJavaWorld	Hello Java World
HelpersMG	Tools for Earth Meteorological Analysis
helsinki	R Tools for Helsinki Open Data
heplots	Visualizing Hypothesis Tests in Multivariate Linear Models
hergm	Hierarchical Exponential-Family Random Graph Models
heritability	Marker-Based Estimation of Heritability Using Individual Plant or Plot Data
hermite	Generalized Hermite Distribution
hett	Heteroscedastic t-regression
het.test	White's Test for Heteroskedasticity
hexbin	Hexagonal Binning Routines
hexView	Viewing Binary Files
hflights	Flights that departed Houston in 2011
hgam	High-dimensional Additive Modelling
hglasso	Learning graphical models with hubs
hglm	Hierarchical Generalized Linear Models
hglm.data	Data for The hglm Package
hgm	Holonomic Gradient Method and Gradient Descent
HGNChelper	Handy Functions for Working With HGNC Gene Symbols and Affymetrix Probeset Identifiers
HH	Statistical Analysis and Data Display: Heiberger and Holland
HHG	Heller-Heller-Gorfine Tests of Independence and Equality of Distributions
hht	The Hilbert-Huang Transform: Tools and Methods
HI	Simulation from distributions supported by nested hyperplanes
HiCfeat	Multiple Logistic Regression for 3D Chromatin Domain Border Analysis
HiClimR	Hierarchical Climate Regionalization
HiCseg	Detection of domains in HiC data
hiddenf	The All-Configurations, Maximum-Interaction F-Test for Hidden Additivity
HiddenMarkov	Hidden Markov Models
HiDimDA	High Dimensional Discriminant Analysis
HiDimMaxStable	Inference on High Dimensional Max-Stable Distributions
hierarchicalDS	Functions For Performing Hierarchical Analysis of Distance Sampling Data
hierband	Convex Banding of the Covariance Matrix
hierDiversity	Hierarchical Multiplicative Partitioning of Complex Phenotypes
hierfstat	Estimation and Tests of Hierarchical F-Statistics
hierNet	A Lasso for Hierarchical Interactions
HierO	A graphical user interface for calculating power and sample size for hierarchical data
hier.part	Hierarchical Partitioning
hiertest	Convex Hierarchical Testing of Interactions
HIest	Hybrid index estimation
highcharter	A Wrapper for the 'Highcharts' Library
highD2pop	Two-Sample Tests for Equality of Means in High Dimension
HighDimOut	Outlier Detection Algorithms for High-Dimensional Data
highfrequency	Tools For Highfrequency Data Analysis
highlight	Syntax Highlighter
highmean	Two-Sample Tests for High-Dimensional Mean Vectors
highr	Syntax Highlighting for R Source Code
highriskzone	Determining and Evaluating High-Risk Zones
highTtest	Simultaneous Critical Values for t-Tests in Very High Dimensions
hillmakeR	Perform occupancy analysis
HiLMM	Estimation of Heritability in Linear Mixed Models
hindexcalculator	H-Index Calculator using Data from a Web of Science (WoS) Citation Report
hint	Tools for hypothesis testing based on Hypergeometric Intersection distributions
HiPLARM	High Performance Linear Algebra in R
hiPOD	hierarchical Pooled Optimal Design
hisemi	Hierarchical Semiparametric Regression of Test Statistics
hisse	Hidden State Speciation and Extinction
HistData	Data Sets from the History of Statistics and Data Visualization
HistDAWass	Histogram-Valued Data Analysis
histmdl	A Most Informative Histogram-Like Model
histogram	Construction of regular and irregular histograms with different options for automatic choice of bins
HistogramTools	Utility Functions for R Histograms
historydata	Data Sets for Historians
hit	Hierarchical Inference Testing
hitandrun	"Hit and Run" and "Shake and Bake" for Sampling Uniformly from Convex Shapes
hive	Hadoop InteractiVE
HiveR	2D and 3D Hive Plots for R
HIV.LifeTables	HIV calibrated model life tables for countries with generalized HIV epidemics
HK80	Conversion Tools for HK80 Geographical Coordinate System
hkevp	Spatial Extreme Value Analysis with the Hierarchical Model of Reich and Shaby (2012)
HKprocess	Hurst-Kolmogorov Process
HLMdiag	Diagnostic Tools for Hierarchical (Multilevel) Linear Models
HLSM	Hierarchical Latent Space Network Model (HLSM)
HMDHFDplus	Read HMD and HFD Data from the Web
hmeasure	The H-measure and other scalar classification performance metrics
Hmisc	Harrell Miscellaneous
HMM	HMM - Hidden Markov Models
HMMCont	Hidden Markov Model for Continuous Observations Processes
hmm.discnp	Hidden Markov models with discrete non-parametric observation distributions
hmmm	hierarchical multinomial marginal models
HMMpa	Analysing accelerometer data using hidden Markov models
HMP	Hypothesis Testing and Power Calculations for Comparing Metagenomic Samples from HMP
HMPTrees	Statistical Object Oriented Data Analysis of RDP-Based Taxonomic Trees from Human Microbiome Data
HMR	Flux Estimation with Static Chamber Data
hnp	Half-Normal Plots with Simulation Envelopes
hoa	Higher Order Likelihood Inference
hoardeR	Information Retrieval for Genetic Datasets
holdem	Texas Holdem simulator
Holidays	Holiday and Half-Day Data, for Use with the 'TimeWarp' Package
homals	Gifi Methods for Optimal Scaling
homeR	Functions useful for building physics
homomorpheR	Homomorphic Computations in R
HomoPolymer	Theoretical Model to Simulate Radical Polymerization
homtest	Homogeneity tests for Regional Frequency Analysis
hopbyhop	Transmissions and Receptions in a Hop by Hop Network
hornpa	Horn's (1965) Test to Determine the Number of Components/Factors
hot.deck	Multiple Hot-Deck Imputation
HotDeckImputation	Hot Deck Imputation Methods for Missing Data
Hotelling	Hotelling's T-squared test and variants
hotspot	Software Hotspot Analysis
hotspots	Hot spots
housingData	U.S. Housing Data from 2008 to 2016
howmany	A lower bound for the number of correct rejections
HPbayes	Heligman Pollard mortality model parameter estimation using Bayesian Melding with Incremental Mixture Importance Sampling
hpcwld	High Performance Cluster Models Based on Kiefer-Wolfowitz Recursion
hpoPlot	Functions for Plotting HPO Terms
hqmisc	Miscellaneous convenience functions and dataset
hqreg	Regularization Paths for Huber Loss Regression and Quantile Regression Penalized by Lasso or Elastic-Net
HRM	High-Dimensional Repeated Measures
hrr	Horizontal rule for the R language
HSAR	Hierarchical Spatial Autoregressive Model (HSAR)
HSAUR	A Handbook of Statistical Analyses Using R (1st Edition)
HSAUR2	A Handbook of Statistical Analyses Using R (2nd Edition)
HSAUR3	A Handbook of Statistical Analyses Using R (3rd Edition)
hsdar	Manage, Analyse and Simulate Hyperspectral Data
hSDM	hierarchical Bayesian species distribution models
hsicCCA	Canonical Correlation Analysis based on Kernel Independence Measures
hsmm	Hidden Semi Markov Models
hsphase	Phasing, Pedigree Reconstruction, Sire Imputation and Recombination Events Identification of Half-sib Families Using SNP Data
HSROC	Meta-Analysis of Diagnostic Test Accuracy when Reference Test is Imperfect
HSSVD	Biclustering with Heterogeneous Variance
htmltab	Assemble Data Frames from HTML Tables
htmlTable	Advanced Tables for Markdown/HTML
htmltools	Tools for HTML
HTMLUtils	Facilitates Automated HTML Report Creation
htmlwidgets	HTML Widgets for R
hts	Hierarchical and Grouped Time Series
HTSCluster	Clustering high throughput sequencing (HTS) data
httk	High-Throughput Toxicokinetics
httpcache	Query Cache for HTTP Clients
httpcode	'HTTP' Status Code Helper
httping	'Ping' 'URLs' to Time 'Requests'
httpRequest	Basic HTTP Request
httpuv	HTTP and WebSocket Server Library
httr	Tools for Working with URLs and HTTP
huge	High-Dimensional Undirected Graph Estimation
HUM	compute HUM value and visualize ROC curves
humanFormat	Human-friendly formatting functions
humaniformat	A Parser for Human Names
humarray	Simplify Analysis and Annotation of Human Microarray Datasets
hunspell	Morphological Analysis and Spell Checker for R
hwde	Models and Tests for Departure from Hardy-Weinberg Equilibrium and Independence Between Loci
HWEBayes	Bayesian investigation of Hardy-Weinberg Equilibrium via estimation and testing
HWEintrinsic	Objective Bayesian Testing for the Hardy-Weinberg Equilibrium Problem
HW.pval	Testing Hardy-Weinberg Equilibrium for Multiallelic Genes
hwriter	HTML Writer - Outputs R objects in HTML format
hwriterPlus	hwriterPlus: Extending the hwriter Package
hwwntest	Tests of White Noise using Wavelets
HWxtest	Exact Tests for Hardy-Weinberg Proportions
hybridEnsemble	Build, Deploy and Evaluate Hybrid Ensembles
hybridHclust	Hybrid Hierarchical Clustering
HybridMC	Implementation of the Hybrid Monte Carlo and Multipoint Hybrid Monte Carlo sampling techniques
hybridModels	Stochastic Hybrid Models in Dynamic Networks
HydeNet	Hybrid Bayesian Networks Using R and JAGS
hydroApps	Tools and models for hydrological applications
hydrogeo	Groundwater data presentation and interpretation
hydroGOF	Goodness-of-fit functions for comparison of simulated and observed hydrological time series
HydroMe	R codes for estimating water retention and infiltration model parameters using experimental data
hydroPSO	Particle Swarm Optimisation, with focus on Environmental Models
hydrostats	Hydrologic Indices for Daily Time Series Data
hydroTSM	Time series management, analysis and interpolation for hydrological modelling
hyfo	Hydrology and Climate Forecasting
HyperbolicDist	The hyperbolic distribution
hyperdirichlet	A generalization of the Dirichlet distribution
hypergea	Hypergeometric tests
hypergeo	The Gauss Hypergeometric Function
hyperSpec	Work with Hyperspectral Data, i.e. Spectra + Meta Information (Spatial, Time, Concentration, ...)
hypervolume	High-Dimensional Kernel Density Estimation and Geometry Operations
hyphenatr	Tools to Hyphenate Strings Using the 'Hunspell' Hyphenation Library
HyPhy	Macroevolutionary phylogentic analysis of species trees and gene trees
hypothesestest	Confidence Intervals and Tests of Statistical Hypotheses
hysteresis	Tools for Modeling Rate-Dependent Hysteretic Processes and Ellipses
hzar	Hybrid Zone Analysis using R
0	
IalsaSynthesis	Synthesizing Information Across Collaborating Research
IASD	Model Selection for Index of Asymmetry Distribution
IAT	Functions to use with data from the Implicit Association Test
IATscores	Implicit Association Test Scores Using Robust Statistics
iBATCGH	Integrative Bayesian Analysis of Transcriptomic and CGH Data
ibd	INCOMPLETE BLOCK DESIGNS
IBDhaploRtools	Functions for the Analysis of IBD Haplo Output
IBDLabels	Convert Between Different IBD-State Labelling Schemes
ibdreg	Regression Methods for IBD Linkage With Covariates
IBDsim	Simulation of Chromosomal Regions Shared by Family Members
ibeemd	Irregular-lattice based ensemble empirical mode decomposition
ibelief	Belief Function Implementation
IBHM	Approximation using the IBHM method
ibmdbR	IBM in-Database Analytics for R
Iboot	Iboot: iterated bootstrap tests and confidence sets
ibr	Iterative Bias Reduction
IBrokers	R API to Interactive Brokers Trader Workstation
iBUGS	An Interface to WinBUGS/OpenBUGS/JAGS by gWidgets
iC10	A Copy Number and Expression-Based Classifier for Breast Tumours
iC10TrainingData	Training datasets for iC10 package
IC2	Inequality and Concentration Indices and Curves
ic50	Standardized high-throughput evaluation of cell-based compound screens
ica	Independent Component Analysis
ICAFF	Imperialist Competitive Algorithm
icamix	Estimation of ICA Mixture Models
icaOcularCorrection	Independent Components Analysis (ICA) based artifact correction
icapca	Mixed ICA/PCA
icarus	Calibrates and Reweights Units in Samples
ICBayes	Bayesian Semiparametric Models for Interval-Censored Data
ICC	Facilitating Estimation of the Intraclass Correlation Coefficient
iccbeta	Multilevel model intraclass correlation for slope heterogeneity
ICC.Sample.Size	Calculation of Sample Size and Power for ICC
icd9	Tools for Working with ICD-9 Codes, and Finding Comorbidities
ICE	Iterated Conditional Expectation
ICEbox	Individual Conditional Expectation Plot Toolbox
ICEinfer	Incremental Cost-Effectiveness (ICE) Statistical Inference from Two Unbiased Samples
icenReg	Regression Models for Interval Censored Data
icensmis	Study Design and Data Analysis in the Presence of Error-Prone Diagnostic Tests and Self-Reported Outcomes
ICGE	Estimation of number of clusters and identification of atypical units
ICGOR	Fit Generalized Odds Rate Hazards Model with Interval Censored Data
ic.infer	Inequality constrained inference in linear normal situations
iClick	A Button-Based GUI for Financial and Economic Data Analysis
iCluster	Integrative clustering of multiple genomic data types
icRSF	A Modified Random Survival Forest Algorithm
ICS	Tools for Exploring Multivariate Data via ICS/ICA
ICSNP	Tools for Multivariate Nonparametrics
ICsurv	A package for semiparametric regression analysis of interval-censored data
icsw	Inverse Compliance Score Weighting
idbg	R debugger
idbr	R Interface to the US Census Bureau International Data Base API
idendr0	Interactive Dendrograms
identifyr	Clean Unique Identifiers
identity	Jacquard Condensed Coefficients of Identity
idm	Incremental Decomposition Methods
IDPmisc	Utilities of Institute of Data Analyses and Process Design (www.idp.zhaw.ch)
IDPSurvival	Imprecise Dirichlet Process for Survival Analysis
idr	Irreproducible discovery rate
IDTurtle	Identify Turtles by their Plastral Biometries
iDynoR	R Analysis package for iDynoMiCS Simulation Results
ieeeround	Functions to set and get the IEEE rounding mode
iemisc	Irucka Embry's Miscellaneous Functions
iemiscdata	Irucka Embry's Miscellaneous Data Collection
ifa	Independent Factor Analysis
iFad	An integrative factor analysis model for drug-pathway association inference
ifaTools	Toolkit for Item Factor Analysis with OpenMx
ifctools	Italian Fiscal Code ('Codice Fiscale') Utilities
IFP	Identifying Functional Polymorphisms
ifs	Iterated Function Systems
ifultools	Insightful Research Tools
iGasso	Statistical Tests and Utilities for Genetic Association
IGM.MEA	IGM MEA Analysis
IgorR	Read Binary Files Saved by 'Igor Pro' (Including 'Neuromatic' Data)
igraph	Network Analysis and Visualization
igraphdata	A Collection of Network Data Sets for the 'igraph' Package
igraphinshiny	Use 'shiny' to Demo 'igraph'
igraphtosonia	Convert iGraph graps to SoNIA .son files
ig.vancouver.2014.topcolour	Instagram 2014 Vancouver Top Colour Dataset
ihs	Inverse Hyperbolic Sine Distribution
iki.dataclim	Consistency, Homogeneity and Summary Statistics of Climatological Data
iLaplace	Improved Laplace Approximation for Integrals of Unimodal Functions
ilc	Lee-Carter Mortality Models using Iterative Fitting Algorithms
IM	Orthogonal Moment Analysis
imager	Image Processing Library Based on CImg
IMak	Item Maker
Imap	Interactive Mapping
imguR	An Imgur.com API Client Package
IMIS	Increamental Mixture Importance Sampling
immer	Item Response Models for Multiple Ratings
IMP	Interactive Model Performance Evaluation
ImpactIV	Identifying Causal Effect for Multi-Component Intervention Using Instrumental Variable Method
imPois	Imprecise Inferential Framework for Poisson Sampling Model
import	An Import Mechanism for R
ImportExport	Import and Export Data
imprProbEst	Minimum distance estimation in an imprecise probability model
imputeLCMD	A collection of methods for left-censored missing data imputation
imputeMDR	The Multifactor Dimensionality Reduction (MDR) Analysis for Incomplete Data
imputeMissings	Impute Missing Values in a Predictive Context
imputeTS	Time Series Missing Value Imputation
imputeYn	Imputing the Last Largest Censored Observation(s) Under Weighted Least Squares
in2extRemes	Into the extRemes Package
inarmix	Mixture models for longitudinal count data
inbreedR	Analysing Inbreeding Based on Genetic Markers
IndependenceTests	Nonparametric tests of independence between random vectors
indicspecies	Relationship Between Species and Groups of Sites
inegiR	Integrate INEGI’s (Mexican Stats Office) API with R
ineq	Measuring Inequality, Concentration, and Poverty
InfDim	Infine-dimensional model (IDM) to analyse phenotypic variation in growth trajectories
inference	Functions to extract inferential values of a fitted model object
InferenceSMR	Inference about the standardized mortality ratio when evaluating the effect of a screening program on survival
inferference	Methods for Causal Inference with Interference
inflection	Finds the inflection point of a curve
influence.ME	Tools for Detecting Influential Data in Mixed Effects Models
influenceR	Software Tools to Quantify Structural Importance of Nodes in a Network
influence.SEM	Case Influence in Structural Equation Models
infoDecompuTE	Information Decomposition of Two-phase Experiments
Information	Data Exploration with Information Theory (Weight-of-Evidence and Information Value)
InformationValue	Performance Analysis and Companion Functions for Binary Classification Models
informR	Sequence Statistics for Relational Event Models
infotheo	Information-Theoretic Measures
infra	An Infrastructure Proxy Function
infuser	A Very Basic Templating Engine
infutil	Information Utility
injectoR	R Dependency Injection
INLABMA	Bayesian Model Averaging with INLA
inline	Functions to Inline C, C++, Fortran Function Calls from R
inlinedocs	Convert inline comments to documentation
inpdfr	Analyse Text Documents Using Ecological Tools
InPosition	Inference Tests for ExPosition
insideRODE	insideRODE includes buildin functions with deSolve solver and C/FORTRAN interfaces to nlme, together with compiled codes
InSilicoVA	Probabilistic Verbal Autopsy Coding with 'InSilicoVA' Algorithm
insol	Solar Radiation
install.load	Check, Install and Load CRAN & USGS GRAN Packages
installr	Using R to Install Stuff (Such As: R, Rtools, RStudio, Git, and More!)
instaR	Access to Instagram API via R
insuranceData	A Collection of Insurance Datasets Useful in Risk Classification in Non-life Insurance
intamap	procedures for automated interpolation
intamapInteractive	procedures for automated interpolation - methods only to be used interactively, not included in intamap package
intcox	Iterated Convex Minorant Algorithm for interval censored event data
IntegratedJM	Joint Modelling of the Gene-Expression and Bioassay Data, Taking Care of the Effect Due to a Fingerprint Feature
Interact	Tests for marginal interactions in a 2 class response model
interactionTest	Calculates Critical Test Statistics to Control False Discovery and Familywise Error Rates in Marginal Effects Plots
InteractiveIGraph	interactive network analysis and visualization
interAdapt	interAdapt
Interatrix	Compute Chi-Square Measures with Corrections
intercure	Cure Rate Estimators for Interval Censored Data
interferenceCI	Exact Confidence Intervals in the Presence of Interference
intergraph	Coercion Routines for Network Data Objects
internetarchive	An API Client for the Internet Archive
interplot	Plot the Effects of Variables in Interaction Terms
Interpol	Interpolation of amino acid sequences
Interpol.T	Hourly interpolation of multiple temperature daily series
interpretR	Binary Classifier and Regression Model Interpretation Functions
InterSIM	Simulation of Inter-Related Genomic Datasets
InterVA4	Replicate and Analyse InterVA4
interval	Weighted Logrank Tests and NPMLE for interval censored data
intervals	Tools for Working with Points and Intervals
interventionalDBN	Interventional Inference for Dynamic Bayesian Networks
IntLik	Numerical Integration for Integrated Likelihood
intpoint	linear programming solver by the interior point method and graphically (two dimensions)
inTrees	Interpret Tree Ensembles
intReg	Interval Regression
intRegGOF	Integrated Regression Goodness of Fit
introgress	methods for analyzing introgression between divergent lineages
intsvy	International Assessment Data Manager
InvariantCausalPrediction	Invariant Causal Prediction
InventorymodelPackage	Inventorymodel
investr	Inverse Estimation/Calibration Functions
invGauss	Threshold regression that fits the (randomized drift) inverse Gaussian distribution to survival data
invLT	Inversion of Laplace-Transformed Functions
io	A Unified Framework for Input-Output Operations in R
ioncopy	Calling Copy Number Alterations in Amplicon Sequencing Data
ionflows	Calculate the Number of Required Flows for Semiconductor Sequencing
ionr	Test for Indifference of Indicator
iosmooth	Functions for smoothing with infinite order flat-top kernels
iotools	I/O Tools for Streaming
ipdmeta	Tools for subgroup analyses with multiple trial data using aggregate statistics
ipdw	Spatial Interpolation by Inverse Path Distance Weighting
ipflasso	Integrative Lasso with Penalty Factors
ipfp	Fast Implementation of the Iterative Proportional Fitting Procedure in C
iplots	iPlots - interactive graphics for R
IPMpack	Builds and analyses Integral Projection Models (IPMs)
ipred	Improved Predictors
ips	Interfaces to Phylogenetic Software in R
IPSUR	Introduction to Probability and Statistics Using R
iptools	Manipulate, Validate and Resolve 'IP' Addresses
ipw	Estimate Inverse Probability Weights
IQCC	Improved Quality Control Charts
iqLearn	Interactive Q-Learning
irace	Iterated Racing Procedures
iRafNet	Integrative Random Forest for Gene Regulatory Network Inference
iRefR	iRefIndex Manager
iRegression	Regression methods for interval-valued variables
iRepro	Reproducibility for Interval-Censored Data
IRISMustangMetrics	Statistics and Metrics for Seismic Data
IRISSeismic	Classes and Methods for Seismic Data Analysis
irlba	Fast Truncated SVD, PCA and Symmetric Eigendecomposition for Large Dense and Sparse Matrices
irr	Various Coefficients of Interrater Reliability and Agreement
irtoys	Simple interface to the estimation and plotting of IRT models
irtProb	Utilities and Probability Distributions Related to Multidimensional Person Item Response Models
irtrees	Estimation of Tree-Based Item Response Models
IRTShiny	Item Response Theory via Shiny
isa2	The Iterative Signature Algorithm
ISBF	Iterative Selection of Blocks of Features - ISBF
isdals	Provides datasets for Introduction to Statistical Data Analysis for the Life Sciences
ISDA.R	interval symbolic data analysis for R
IsingFit	Fitting Ising models using the eLasso method
isingLenzMC	Monte Carlo for Classical Ising Model
IsingSampler	Sampling Methods and Distribution Functions for the Ising Model
ISLR	Data for An Introduction to Statistical Learning with Applications in R
ismev	An Introduction to Statistical Modeling of Extreme Values
Iso	Functions to Perform Isotonic Regression
IsoCI	Confidence intervals for current status data based on transformations and bootstrap
isocir	Isotonic Inference for Circular data
ISOcodes	Selected ISO Codes
IsoGene	Order-Restricted Inference for Microarray Experiments
isopam	Isopam (Clustering)
isopat	Calculation of isotopic pattern for a given molecular formula
isoph	Isotonic Proportional Hazards Model
ISOpureR	Deconvolution of Tumour Profiles
isotone	Active Set and Generalized PAVA for Isotone Optimization
isotonic.pen	Penalized Isotonic Regression in one and two dimensions
IsotopeR	Stable Isotope Mixing Model
ISOweek	Week of the year and weekday according to ISO 8601
isva	Independent Surrogate Variable Analysis
ISwR	Introductory Statistics with R
itcSegment	Individual Tree Crowns Segmentation
ITEMAN	Classical Item Analysis
iteRates	Parametric rate comparison
iterators	Provides Iterator Construct for R
iterLap	Approximate probability densities by iterated Laplace Approximations
iterpc	Efficient Iterator for Permutations and Combinations
itertools	Iterator Tools
itertools2	itertools2: Functions creating iterators for efficient looping
itree	Tools for classification and regression trees, with an emphasis on interpretability
its	Irregular Time Series
itsadug	Interpreting Time Series and Autocorrelated Data Using GAMMs
itsmr	Time series analysis package for students
IUPS	Incorporating Uncertainties in Propensity Scores
ivbma	Bayesian Instrumental Variable Estimation and Model Determination via Conditional Bayes Factors
ivfixed	Instrumental fixed effect panel data model
ivlewbel	Uses heteroscedasticity to estimate mismeasured and endogenous regressor models
ivmodel	Statistical Inference and Sensitivity Analysis for Instrumental Variables Model
ivpack	Instrumental Variable Estimation
ivpanel	Instrumental Panel Data Models
ivprobit	Instrumental variables probit model
iWeigReg	Improved methods for causal inference and missing data problems
iWISA	Wavelet-Based Index of Storm Activity
0	
jaatha	Simulation-Based Maximum Likelihood Parameter Estimation
jackknifeKME	Jackknife Estimates of Kaplan-Meier Estimators or Integrals
jackstraw	Statistical Inference of Variables Driving Systematic Variation
JacobiEigen	Classical Jacobi Eigensolution Algorithm
JADE	Blind Source Separation Methods Based on Joint Diagonalization and Some BSS Performance Criteria
jagsUI	A Wrapper Around 'rjags' to Streamline 'JAGS' Analyses
JAGUAR	Joint Analysis of Genotype and Group-Specific Variability Using a Novel Score Test Approach to Map Expression Quantitative Trait Loci (eQTL)
james.analysis	Analysis Tools for the 'JAMES' Framework
JASPAR	R modules for JASPAR databases: a collection of transcription factor DNA-binding preferences, modeled as matrices
JavaGD	Java Graphics Device
JBTools	Misc Small Tools and Helper Functions for Other Code of J. Buttlar
jetset	One-to-One Gene-Probeset Mapping for Affymetrix Human Microarrays
JGEE	Joint Generalized Estimating Equation Solver
JGL	Performs the Joint Graphical Lasso for sparse inverse covariance estimation on multiple classes
JGR	JGR - Java GUI for R
jiebaR	Chinese Text Segmentation
jiebaRD	Chinese Text Segmentation Data for jiebaR Package
JM	Joint Modeling of Longitudinal and Survival Data
JMbayes	Joint Modeling of Longitudinal and Time-to-Event Data under a Bayesian Approach
jmcm	Joint Mean-Covariance Models using 'Armadillo' and S4
JMdesign	Joint Modeling of Longitudinal and Survival Data - Power Calculation
jmetrik	Tools for Interacting with 'jMetrik'
Jmisc	Julian Miscellaneous Function
jmotif	Tools for Time Series Analysis Based on Symbolic Aggregate Dicretization
jmuOutlier	Permutation Tests for Nonparametric Statistics
Johnson	Johnson Transformation
JohnsonDistribution	Johnson Distribution
joineR	Joint modelling of repeated measurements and time-to-event data
joint.Cox	Penalized Likelihood Estimation and Prediction under the Joint Cox Models Between Tumour Progression and Death for Meta-Analysis
jointDiag	Joint Approximate Diagonalization of a set of square matrices
jointPm	Risk estimation using the joint probability method
JointRegBC	Joint Modelling of Mixed Correlated Binary and Continuous Responses : A Latent Variable Approach
jomo	Multilevel Joint Modelling Multiple Imputation
JOP	Joint Optimization Plot
JoSAE	Functions for some Unit-Level Small Area Estimators and their Variances
jpeg	Read and write JPEG images
JPEN	Covariance and Inverse Covariance Matrix Estimation Using Joint Penalty
JPSurv	Methods for population-based cancer survival analysis
jqr	Client for 'jq', a JSON Processor
JRF	Joint Random Forest (JRF) for the Simultaneous Estimation of Multiple Related Networks
jrich	Jack-Knife Support for Evolutionary Distinctiveness Indices I and W
jrvFinance	Basic Finance; NPV/IRR/Annuities/Bond-Pricing; Black Scholes
js	Tools for Working with JavaScript in R
jSonarR	jSonar Analytics Platform API for R
jsonlite	A Robust, High Performance JSON Parser and Generator for R
jtrans	Johnson Transformation for Normality
jug	Create a Simple Web API for your R Functions
Julia	Fractal Image Data Generator
junr	Access Open Data Through the Junar API
jvnVaR	Value at Risk
0	
KANT	Package to identify and sort genes overexpressed
kantorovich	Kantorovich Distance Between Probability Measures
KappaGUI	GUI for Cohen's and Fleiss' Kappa
kappalab	Non-Additive Measure and Integral Manipulation Functions
kappaSize	Sample Size Estimation Functions for Studies of Interobserver Agreement
KappaV	Calculates "vectorial Kappa", an index of congruence between patchy mosaics
kaps	K-Adaptive Partitioning for Survival data
karaoke	Remove Vocals from a Song
KATforDCEMRI	Kinetic analysis and visualization of DCE-MRI data
kcirt	k-Cube Thurstonian IRT Models
kdecopula	Kernel Smoothing for Bivariate Copula Densities
kdetrees	Nonparametric method for identifying discordant phylogenetic trees
kedd	Kernel Estimator and Bandwidth Selection for Density and Its Derivatives
keep	Arrays with Better Control over Dimension Dropping
kelvin	Calculate Solutions to the Kelvin Differential Equation using Bessel Functions
Kendall	Kendall rank correlation and Mann-Kendall trend test
kequate	The Kernel Method of Test Equating
kerdiest	Nonparametric kernel estimation of the distribution function. Bandwidth selection and estimation of related functions
KERE	Expectile Regression in Reproducing Kernel Hilbert Space
kergp	Gaussian Process Laboratory
kernDeepStackNet	Kernel Deep Stacking Networks
kerndwd	Distance Weighted Discrimination (DWD) and Kernel Methods
kernelFactory	Kernel Factory: An Ensemble of Kernel Machines
Kernelheaping	Kernel Density Estimation for Heaped and Rounded Data
kernlab	Kernel-Based Machine Learning Lab
KernSmooth	Functions for Kernel Smoothing Supporting Wand & Jones (1995)
KernSmoothIRT	Nonparametric Item Response Theory
keyplayer	Locating Key Players in Social Networks
keypress	Wait for a Key Press in a Terminal
KFAS	Kalman Filter and Smoother for Exponential Family State Space Models
kfigr	Integrated Code Chunk Anchoring and Referencing for R Markdown Documents
KFKSDS	Kalman Filter, Smoother and Disturbance Smoother
kimisc	Kirill's Miscellaneous Functions
kin.cohort	Analysis of Kin-Cohort Studies
kineticF	Framework for the Analysis of Kinetic Visual Field Data
kinfit	Routines for Fitting Kinetic Models to Chemical Degradation Data
kinship2	Pedigree Functions
kissmig	a Keep It Simple Species Migration Model
kitagawa	Spectral response of water wells to harmonic strain and pressure
kknn	Weighted k-Nearest Neighbors
klaR	Classification and visualization
klausuR	Multiple Choice Test Evaluation
klin	Linear equations with Kronecker structure
kmc	Kaplan-Meier Estimator with Constraints for Right Censored Data – a Recursive Computational Algorithm
km.ci	Confidence intervals for the Kaplan-Meier estimator
kmconfband	Kaplan-Meier Simultaneous Confidence Band for the Survivor Function
KMDA	Kernel-Based Metabolite Differential Analysis
kmeans.ddR	Distributed k-Means for Big Data using 'ddR' API
kmi	Kaplan-Meier multiple imputation for the analysis of cumulative incidence functions in the competing risks setting
Kmisc	Kevin Miscellaneous
kml	K-Means for Longitudinal Data
kml3d	K-Means for Joint Longitudinal Data
kmlcov	Clustering longitudinal data using the likelihood as a metric of distance
kmlShape	K-Means for Longitudinal Data using Shape-Respecting Distance
kmodR	K-Means with Simultaneous Outlier Detection
KMsurv	Data sets from Klein and Moeschberger (1997), Survival Analysis
knitcitations	Citations for 'Knitr' Markdown Files
knitLatex	'Knitr' Helpers - Mostly Tables
knitr	A General-Purpose Package for Dynamic Report Generation in R
knitrBootstrap	Knitr Bootstrap Framework
knncat	Nearest-neighbor Classification with Categorical Variables
knnGarden	Multi-distance based k-Nearest Neighbors
knnIndep	Independence tests and benchmarks
knockoff	Knockoff Filter for Controlling the False Discovery Rate
kobe	Tools for the provision of scientific fisheries management advice
KODAMA	Knowledge discovery by accuracy maximization
kofnGA	A Genetic Algorithm for Fixed-Size Subset Selection
KOGMWU	Functional Summary and Meta-Analysis of Gene Expression Data
kohonen	Supervised and Unsupervised Self-Organising Maps
kolmim	An Improved Evaluation of Kolmogorov's Distribution
KoNLP	Korean NLP Package
koRpus	An R Package for Text Analysis
KoulMde	Koul's Minimum Distance Estimation in Linear Regression and Autoregression Model
Kpart	Spline Fitting
kpodclustr	Method for Clustering Partially Observed Data
kriens	Continuation Passing Style Development
kriging	Ordinary Kriging
KrigInv	Kriging-based Inversion for Deterministic and Noisy Computer Experiments
KRLS	Kernel-based Regularized Least squares (KRLS)
krm	Kernel Based Regression Models
ks	Kernel Smoothing
kSamples	K-Sample Rank Tests and their Combinations
KScorrect	Lilliefors-Corrected Kolmogorov-Smirnoff Goodness-of-Fit Tests
kselection	Selection of K in K-Means Clustering
ksrlive	Identify Kinase Substrate Relationships Using Dynamic Data
kst	Knowledge Space Theory
ktsolve	Configurable function for solving families of nonlinear equations
ktspair	k-Top Scoring Pairs for Microarray Classification
kulife	Datasets and functions from the (now non-existing) Faculty of Life Sciences, University of Copenhagen
kwb.hantush	Calculation of Groundwater Mounding Beneath an Infiltration Basin
kyotil	Utility Functions by Youyi, Krisz and Others
kza	Kolmogorov-Zurbenko Adaptive Filters
kzft	Kolmogorov-Zurbenko Fourier Transform and Applications
kzs	Kolmogorov-Zurbenko Spatial Smoothing and Applications
0	
L1pack	Routines for L1 Estimation
l2boost	l2boost - Friedman's boosting algorithm for regularized linear regression
labdsv	Ordination and Multivariate Analysis for Ecology
labeledLoop	Labeled Loop
labeling	Axis Labeling
labelrank	Predicting Rankings of Labels
label.switching	Relabelling MCMC Outputs of Mixture Models
labeltodendro	Convert labels or tables to a dendrogram
labstatR	Libreria Del Laboratorio Di Statistica Con R
labstats	Data Sets for the Book "Experimental Design for Laboratory Biologists"
laeken	Estimation of indicators on social exclusion and poverty
laercio	Duncan test, Tukey test and Scott-Knott test
LaF	Fast access to large ASCII files
laGP	Local Approximate Gaussian Process Regression
Lahman	Sean Lahman's Baseball Database
LakeMetabolizer	Tools for the Analysis of Ecosystem Metabolism
lakemorpho	Lake morphometry in R
laketemps	Lake Temperatures Collected by Situ and Satellite Methods from 1985-2009
Lambda4	Collection of Internal Consistency Reliability Coefficients
lambda.r	Modeling Data with Functional Programming
lambda.tools	Tools for Modeling Data with Functional Programming
LambertW	Probabilistic Models to Analyze and Gaussianize Heavy-Tailed, Skewed Data
lamW	Lambert-W Function
LANDD	Liquid Association for Network Dynamics Detection
landest	Landmark Estimation of Survival and Treatment Effect
landpred	Landmark Prediction of a Survival Outcome
landsat	Radiometric and topographic correction of satellite imagery
landsat8	Landsat 8 Imagery Rescaled to Reflectance, Radiance and/or Temperature
Langevin	Langevin Analysis in One and Two Dimensions
languageR	Data sets and functions with "Analyzing Linguistic Data: A practical introduction to statistics"
LaplaceDeconv	Laplace Deconvolution with Noisy Discrete Non-Equally Spaced Observations on a Finite Time Interval
lar	History of labour relations package
LARF	Local Average Response Functions for Instrumental Variable Estimation of Treatment Effects
lars	Least Angle Regression, Lasso and Forward Stagewise
laser	Likelihood Analysis of Speciation/Extinction Rates from Phylogenies
lasso2	L1 constrained estimation aka ‘lasso’
lassoscore	High-Dimensional Inference with the Penalized Score Test
lassoshooting	L1 regularized regression (Lasso) solver using the Cyclic Coordinate Descent algorithm aka Lasso Shooting
lasvmR	A Simple Wrapper for the LASVM Solver
latdiag	Draws Diagrams Useful for Checking Latent Scales
latentnet	Latent Position and Cluster Models for Statistical Networks
Laterality	Functions to Calculate Common Laterality Statistics in Primatology
latex2exp	Use LaTeX Expressions in Plots
lattice	Trellis Graphics for R
latticeDensity	Density estimation and nonparametric regression on irregular regions
latticeExtra	Extra Graphical Utilities Based on Lattice
LatticeKrig	Multiresolution Kriging Based on Markov Random Fields
lava	Latent Variable Models
lavaan	Latent Variable Analysis
lavaan.shiny	Latent Variable Analysis with Shiny
lavaan.survey	Complex Survey Structural Equation Modeling (SEM)
lava.tobit	LVM with censored and binary outcomes
lawn	Client for 'Turfjs' for 'Geospatial' Analysis
lawstat	Tools for Biostatistics, Public Policy, and Law
lazy	Lazy Learning for Local Regression
lazyData	A LazyData Facility
lazyeval	Lazy (Non-Standard) Evaluation
lazysql	Lazy SQL Programming
lazyWeave	LaTeX Wrappers for R Users
lba	Latent Budget Analysis for Compositional Data
lbfgs	Limited-memory BFGS Optimization
lbfgsb3	Limited Memory BFGS Minimizer with Bounds on Parameters
lbiassurv	Length-biased correction to survival curve estimation
LCA	Localised Co-Dependency Analysis
LCAextend	Latent Class Analysis (LCA) with familial dependence in extended pedigrees
lcda	Latent Class Discriminant Analysis
LCFdata	Data sets for package “LMERConvenienceFunctions”
LCMCR	Bayesian Nonparametric Latent-Class Capture-Recapture
lcmm	Extended Mixed Models Using Latent Classes and Latent Processes
lcopula	Liouville Copulas
lctools	Local Correlation, Spatial Inequalities, Geographically Weighted Regression and Other Tools
lda	Collapsed Gibbs Sampling Methods for Topic Models
ldamatch	Multivariate Condition Matching by Backwards Elimination Using Linear Discriminant Analysis
ldatuning	Tuning of the LDA Models Parameters
LDAvis	Interactive Visualization of Topic Models
ldbounds	Lan-DeMets Method for Group Sequential Boundaries
LDcorSV	Linkage disequilibrium corrected by the structure and the relatedness
LDheatmap	Graphical display of pairwise linkage disequilibria between SNPs
ldlasso	LD LASSO Regression for SNP Association Study
LDOD	Finding Locally D-optimal optimal designs for some nonlinear and generalized linear models
LDPD	Probability of Default Calibration
ldr	Methods for likelihood-based dimension reduction in regression
LDRTools	Tools for Linear Dimension Reduction
LDtests	Exact tests for Linkage Disequilibrium and Hardy-Weinberg Equilibrium
leaderCluster	Leader Clustering Algorithm
LeafAngle	Analysis and Visualization of Plant Leaf Angle Distributions
LeafArea	Rapid Digital Image Analysis of Leaf Area
leaflet	Create Interactive Web Maps with the JavaScript 'Leaflet' Library
leafletR	Interactive Web-Maps Based on the Leaflet JavaScript Library
LEAPFrOG	Likelihood Estimation of Admixture in Parents From Offspring Genotypes
leapp	latent effect adjustment after primary projection
leaps	regression subset selection
LearnBayes	Functions for Learning Bayesian Inference
learningr	Data and functions to accompany the book "Learning R"
learNN	Examples of Neural Networks
learnstats	An Interactive Environment for Learning Statistics
lefse	Phylogenetic and Functional Analyses for Ecology
leiv	Bivariate Linear Errors-In-Variables Estimation
LeLogicielR	Functions and datasets to accompany the book "Le logiciel R: Maitriser le langage, Effectuer des analyses statistiques" (french)
lessR	Less Code, More Results
lestat	A package for LEarning STATistics
letsR	Tools for Data Handling and Analysis in Macroecology
lettercase	Utilities for Formatting Strings with Consistent Capitalization, Word Breaks and White Space
LexisPlotR	Plot Lexis Diagrams for Demographic Purposes
lfactors	Factors with Levels
lfda	Local Fisher Discriminant Analysis
LFDR.MLE	Estimation of the Local False Discovery Rates by Type II Maximum Likelihood Estimation
lfe	Linear Group Fixed Effects
lfl	Linguistic Fuzzy Logic
lfstat	Calculation of Low Flow Statistics for Daily Stream Flow Data
lga	Tools for linear grouping analysis (LGA)
lgarch	Simulation and Estimation of Log-GARCH Models
lgcp	Log-Gaussian Cox Process
LGEWIS	Tests for Genetic Association/Gene-Environment Interaction in Longitudinal Gene-Environment-Wide Interaction Studies
LGRF	Set-Based Tests for Genetic Association in Longitudinal Studies
lgtdl	A set of methods for longitudinal data objects
lhs	Latin Hypercube Samples
libamtrack	Computational Routines for Proton and Ion Radiotherapy
LiblineaR	Linear Predictive Models Based on the LIBLINEAR C/C++ Library
LiblineaR.ACF	Linear Classification with Online Adaptation of Coordinate Frequencies
Libra	Linearized Bregman Algorithms for Generalized Linear Models
LICORS	Light Cone Reconstruction of States - Predictive State Estimation From Spatio-Temporal Data
LICurvature	Sensitivity Analysis for Case Weight in Normal Linear Regression
lifecontingencies	Financial and Actuarial Mathematics for Life Contingencies
LifeHist	Life History Models of Individuals
LifeTables	Two-Parameter HMD Model Life Table System
lift	Compute the Top Decile Lift and Plot the Lift Curve
liftr	Dockerize R Markdown Documents
LightningR	Tools for Communication with Lightning-Viz Server
lightsout	Implementation of the 'Lights Out' Puzzle Game
LIHNPSD	Poisson Subordinated Distribution
likelihood	Methods for Maximum Likelihood Estimation
likelihoodAsy	Functions for Likelihood Asymptotics
likeLTD	Tools to Evaluate DNA Profile Evidence
likert	Functions to Analyze and Visualize Likert Type Items
LIM	Linear Inverse Model examples and solution methods
limitplot	Jitter/CI Plot with Ordered Points Below the Limit of Detection
limSolve	Solving Linear Inverse Models
linbin	Binning and Plotting of Linearly Referenced Data
LinCal	Static Univariate Frequentist and Bayesian Linear Calibration
LindenmayeR	Functions to Explore L-Systems (Lindenmayer Systems)
LinearizedSVR	Linearized Support Vector Regression
LinearRegressionMDE	Minimum Distance Estimation in Linear Regression Model
linERR	Linear Excess Relative Risk Model
lineup	Lining Up Two Sets of Measurements
linkcomm	Tools for Generating, Visualizing, and Analysing Link Communities in Networks
LinkedMatrix	Column-Linked and Row-Linked Matrices
linkim	Linkage information based genotype imputation method
linkR	3D Lever and Linkage Mechanism Modeling
linLIR	linear Likelihood-based Imprecise Regression
linprog	Linear Programming / Optimization
LinRegInteractive	Interactive Interpretation of Linear Regression Models
LINselect	Selection of Linear Estimators
lint	Tools to check R code style
lintools	Manipulation of Linear Systems of (in)Equalities
lintr	Static R Code Analysis
lira	LInear Regression in Astronomy
liso	Fitting lasso penalised additive isotone models
lisp	List-processing à la SRFI-1
lisrelToR	Import output from LISREL into R
list	Statistical Methods for the Item Count Technique and List Experiment
listenv	Environments Behaving (Almost) as Lists
LIStest	Tests of independence based on the Longest Increasing Subsequence
listWithDefaults	List with Defaults
littler	R at the Command-Line via 'r'
livechatR	R Wrapper for LiveChat REST API
llama	Leveraging Learning to Automatically Manage Algorithms
lle	Locally linear embedding
lllcrc	Local Log-linear Models for Capture-Recapture
LLSR	Data Analysis of Liquid-Liquid Systems
lm.beta	Add Standardized Regression Coefficients to lm-Objects
lm.br	Linear Model with Breakpoint
lme4	Linear Mixed-Effects Models using 'Eigen' and S4
lmec	Linear Mixed-Effects Models with Censored Responses
lmeNB	Compute the Personalized Activity Index Based on a Negative Binomial Model
lmeNBBayes	Compute the Personalized Activity Index Based on a Flexible Bayesian Negative Binomial Model
lmenssp	Linear Mixed Effects Models with Non-Stationary Stochastic Processes
LMERConvenienceFunctions	Model Selection and Post-hoc Analysis for (G)LMER Models
lmerTest	Tests in Linear Mixed Effects Models
lmeSplines	Add smoothing spline modelling capability to nlme
LMest	Latent Markov Models with and without Covariates
lmeVarComp	Testing for a subset of variance components in linear mixed models
lmf	Functions for estimation and inference of selection in age-structured populations
lmfor	Functions for Forest Biometrics
lmm	Linear Mixed Models
lmmlasso	Linear mixed-effects models with Lasso
lmmot	Multiple Ordinal Tobit (MOT) Model
lmms	Linear Mixed Effect Model Splines for Modelling and Analysis of Time Course Data
lmodel2	Model II Regression
lmom	L-moments
lmomco	L-Moments, Censored L-Moments, Trimmed L-Moments, L-Comoments, and Many Distributions
Lmoments	L-Moments and Quantile Mixtures
lmomRFA	Regional frequency analysis using L-moments
lmSupport	Support for Linear Models
lmtest	Testing Linear Regression Models
LncMod	Predicting Modulator and Functional/Survival Analysis
loa	Lattice Options and Add-Ins
localdepth	Local Depth
localgauss	Estimating Local Gaussian Parameters
localsolver	R API to LocalSolver
locfdr	Computes Local False Discovery Rates
LocFDRPois	Functions for Performing Local FDR Estimation when Null and Alternative are Poisson
locfit	Local Regression, Likelihood and Density Estimation
locits	Test of stationarity and localized autocovariance
Lock5Data	Datasets for "Statistics: UnLocking the Power of Data"
Lock5withR	Datasets for 'Statistics: Unlocking the Power of Data'
locpol	Kernel local polynomial regression
lodGWAS	Genome-Wide Association Analysis of a Biomarker Accounting for Limit of Detection
loe	Local Ordinal Embedding
log4r	A simple logging system for R, based on log4j
logbin	Relative Risk Regression Using the Log-Binomial Model
LogConcDEAD	Log-concave Density Estimation in Arbitrary Dimensions
logconcens	Maximum likelihood estimation of a log-concave density based on censored data
logcondens	Estimate a Log-Concave Probability Density from iid Observations
logcondens.mode	Compute MLE of Log-Concave Density on R with Fixed Mode, and Perform Inference for the Mode
logcondiscr	Estimate a Log-Concave Probability Mass Function from Discrete i.i.d. Observations
logconPH	CoxPH Model with Log Concave Baseline Distribution
logging	R logging package
LogicForest	Logic Forest
LOGICOIL	LOGICOIL: multi-state prediction of coiled-coil oligomeric state
LogicReg	Logic Regression
logistf	Firth's bias reduced logistic regression
LogisticDx	Diagnostic Tests for Models with a Binomial Response
logisticPCA	Binary Dimensionality Reduction
LOGIT	Functions, Data and Code for Binary and Binomial Data
logitchoice	Fitting l2-regularized logit choice models via generalized gradient descent
LogitNet	Infer network based on binary arrays using regularized logistic regression
logitnorm	Functions for the logitnormal distribution
loglognorm	Double log normal distribution functions
logmult	Log-Multiplicative Models, Including Association Models
LogrankA	Logrank Test for Aggregated Survival Data
logspline	Logspline Density Estimation Routines
lokern	Kernel Regression Smoothing with Local or Global Plug-in Bandwidth
lomb	Lomb-Scargle Periodogram
longCatEDA	Package for Plotting Categorical Longitudinal and Time-Series Data
longclust	Model-Based Clustering and Classification for Longitudinal Data
longitudinal	Analysis of Multiple Time Course Data
longitudinalData	Longitudinal Data
longmemo	Statistics for Long-Memory Processes (Jan Beran) – Data and Functions
longpower	Sample size calculations for longitudinal data
longurl	Expand Short URLs Using the 'LongURL' API
loo	Efficient Leave-One-Out Cross-Validation and WAIC for Bayesian Models
lookupTable	Look-Up Tables using S4
loop	loop decomposition of weighted directed graphs for life cycle analysis, providing flexbile network plotting methods, and analyzing food chain properties in ecology
LoopAnalyst	A Collection of Tools to Conduct Levins' Loop Analysis
loopr	Uses an Archive to Amend Previous Stages of a Pipe using Current Output
lordif	Logistic Ordinal Regression Differential Item Functioning using IRT
lorec	LOw Rand and sparsE Covariance matrix estimation
LOST	Missing morphometric data simulation and estimation
LotkasLaw	Runs Lotka's Law which is One of the Special Applications of Zipf's Law
LowRankQP	Low Rank Quadratic Programming
lpbrim	LP-BRIM Bipartite Modularity
lpc	Lassoed principal components for testing significance of features
LPCM	Local Principal Curve Methods
lpint	Local polynomial estimators of intensity function or its derivatives
LPM	Linear Parametric Models Applied to Hydrological Series
lpme	Local Polynomial Estimator in Measurement Error Models
LPmerge	Merging linkage maps by linear programming
lpmodeler	Modeler for linear programs (LP) and mixed integer linear programs (MILP)
LPR	Lasso and Partial Ridge
lpridge	Local Polynomial (Ridge) Regression
LPS	Linear Predictor Score, for Binary Inference from Multiple Continuous Variables
lpSolve	Interface to 'Lp_solve' v. 5.5 to Solve Linear/Integer Programs
lpSolveAPI	R Interface to 'lp_solve' Version 5.5.2.0
LPStimeSeries	Learned Pattern Similarity and Representation for Time Series
LPTime	LP Nonparametric Approach to Non-Gaussian Non-Linear Time Series Modelling
lqa	Penalized Likelihood Inference for GLMs
lqmm	Linear Quantile Mixed Models
lqr	Robust Linear Quantile Regression
LRcontrast	Dose Response Signal Detection under Model Uncertainty
lrequire	Sources an R "Module" with Caching & Encapsulation, Returning Exported Vars
lrgs	Linear Regression by Gibbs Sampling
lrmest	Different types of estimators to deal with multicollinearity
LRTH	A Likelihood Ratio Test Accounting for Genetic Heterogeneity
LS2W	Locally stationary two-dimensional wavelet process estimation scheme
LS2Wstat	A Multiscale Test of Spatial Stationarity for LS2W processes
lsa	Latent Semantic Analysis
LSAfun	Applied Latent Semantic Analysis (LSA) Functions
lsbclust	Least-Squares Bilinear Clustering for Three-Way Data
LSC	Local Statistical Complexity - Automatic Pattern Discovery in Spatio-Temporal Data
LSD	Lots of Superior Depictions
LSDinterface	Reading LSD Results (.res) Files
lsdv	Least square dummy variable regression
lsei	Solving Least Squares Problems under Equality/Inequality Constraints
lsgl	Linear Sparse Group Lasso
lshorth	The Length of the Shorth
lsl	Latent Structure Learning
lsmeans	Least-Squares Means
LSMonteCarlo	American options pricing with Least Squares Monte Carlo method
lspls	LS-PLS Models
lsr	Companion to "Learning Statistics with R"
lss	the accelerated failure time model to right censored data based on least-squares principle
LSTS	Locally Stationary Time Series
ltbayes	Simulation-Based Bayesian Inference for Latent Traits of Item Response Models
ltm	Latent Trait Models under IRT
ltmle	Longitudinal Targeted Maximum Likelihood Estimation
LTPDvar	LTPD and AOQL Plans for Acceptance Sampling Inspection by Variables
LTR	Perform LTR analysis on microarray data
ltsa	Linear Time Series Analysis
ltsbase	Ridge and Liu Estimates based on LTS (Least Trimmed Squares) Method
ltsk	Local Time Space Kriging
lubridate	Make Dealing with Dates a Little Easier
luca	Likelihood inference from case-control data Under Covariate Assumptions (LUCA)
lucid	Printing Floating Point Numbers in a Human-Friendly Format
lucr	Currency Formatting and Conversion
lulcc	Land Use Change Modelling in R
lumendb	Lumen Database API Client
Luminescence	Comprehensive Luminescence Dating Data Analysis
lunar	Lunar Phase & Distance, Seasons and Other Environmental Factors
luzlogr	Lightweight Logging for R Scripts
lvm4net	Latent Variable Models for Networks
LVMMCOR	A Latent Variable Model for Mixed Continuous and Ordinal Responses
LW1949	An Automated Approach to Evaluating Dose-Effect Experiments Following Litchfield and Wilcoxon (1949)
lxb	Fast LXB File Reader
lymphclon	Accurate Estimation of Clonal Coincidences and Abundances from Biological Replicates
0	
M3	Reading M3 files
M4comp	Data from the M4 Time Series Forecasting Competition
m4fe	Models for Financial Economics
maboost	Binary and Multiclass Boosting Algorithms
MAc	Meta-Analysis with Correlations
MAclinical	Class prediction based on microarray data and clinical parameters
MAd	Meta-Analysis with Mean Differences
mada	Meta-Analysis of Diagnostic Accuracy
maddison	Maddison Project Database
madness	Automatic Differentiation of Multivariate Operations
mads	Multi-Analysis Distance Sampling
madsim	A Flexible Microarray Data Simulation Model
Maeswrap	Wrapper Functions for MAESTRA/MAESPA
magclass	Data Class and Tools for Handling Spatial-Temporal Data
magic	create and investigate magic squares
magicaxis	Pretty Scientific Plotting with Minor-Tick and log Minor-Tick Support
magrittr	A Forward-Pipe Operator for R
maGUI	A Graphical User Interface for Microarray Data Analysis
mail	Sending Email Notifications from R
mailR	A Utility to Send Emails from R
MAINT.Data	Model and Analyse Interval Data
MakefileR	Create 'Makefiles' Using R
makeProject	Creates an empty package framework for the LCFD format
MALDIquant	Quantitative Analysis of Mass Spectrometry Data
MALDIquantForeign	Import/Export Routines for MALDIquant
mallet	A wrapper around the Java machine learning tool MALLET
MAMA	Meta-Analysis of MicroArray
MAMS	Designing Multi-Arm Multi-Stage Studies
MAMSE	Calculation of Minimum Averaged Mean Squared Error (MAMSE) weights
managelocalrepo	Manage a CRAN-Style Local Repository
MANCIE	Matrix Analysis and Normalization by Concordant Information Enhancement
mangoTraining	Mango Solutions Training Datasets
manifestoR	Access and Process Data and Documents of the Manifesto Project
manipulate	Interactive Plots for RStudio
ManlyMix	Manly Mixture Modeling and Model-Based Clustering
ManyTests	Multiple testing procedures of Cox (2011) and Wang and Cox (2007)
Map2NCBI	Mapping Markers to the Nearest Genomic Feature
MAPA	Multiple Aggregation Prediction Algorithm
mapdata	Extra Map Databases
mapfit	A Tool for PH/MAP Parameter Estimation
MapGAM	Mapping Smoothed Effect Estimates from Individual-Level Data
MAPLES	Smoothed age profile estimation
mapmisc	Utilities for Producing Maps
mapplots	Data Visualisation on Maps
mapproj	Map Projections
mapr	'Visualize' Species Occurrence Data
maps	Draw Geographical Maps
mapStats	Geographic Display of Survey Data Statistics
maptools	Tools for Reading and Handling Spatial Objects
maptpx	MAP Estimation of Topic Models
maptree	Mapping, pruning, and graphing tree models
mapview	Interactive Viewing of Spatial Objects in R
mAr	Multivariate AutoRegressive analysis
42430	Multivariate Autoregressive Modeling for Analysis of Community Time-Series Data
mar1s	Multiplicative AR(1) with Seasonal Processes
marelac	Tools for Aquatic Sciences
MareyMap	Estimation of Meiotic Recombination Rates Using Marey Maps
marg	Approximate marginal inference for regression-scale models
markdown	'Markdown' Rendering for R
marked	Mark-Recapture Analysis for Survival and Abundance Estimation
maRketSim	Market simulator for R
markmyassignment	Automatic Marking of R Assignments
markophylo	Markov Chain Models for Phylogenetic Trees
markovchain	Easy Handling Discrete Time Markov Chains
MarkowitzR	Statistical Significance of the Markowitz Portfolio
marl	Multivariate Analysis Based on Relative Likelihoods
marmap	Import, Plot and Analyze Bathymetric and Topographic Data
marqLevAlg	An algorithm for least-squares curve fitting
MARSS	Multivariate Autoregressive State-Space Modeling
maSAE	Mandallaz' Model-Assisted Small Area Estimators
MASS	Support Functions and Datasets for Venables and Ripley's MASS
MASSTIMATE	Body Mass Estimation Equations for Vertebrates
MasterBayes	ML and MCMC Methods for Pedigree Reconstruction and Analysis
MAT	Multidimensional Adaptive Testing
MATA	Model-Averaged Tail Area Wald (MATA-Wald) Confidence Interval
Matching	Multivariate and Propensity Score Matching with Balance Optimization
MatchingFrontier	Computation of the Balance - Sample Size Frontier in Matching Methods for Causal Inference
matchingMarkets	Analysis of Stable Matchings
matchingR	Matching Algorithms in R and C++
MatchIt	MatchIt: Nonparametric Preprocessing for Parametric Casual Inference
MatchLinReg	Combining Matching and Linear Regression for Causal Inference
matchMulti	Optimal Multilevel Matching using a Network Algorithm
matconv	A Code Converter from the Matlab/Octave Language to R
mateable	Tools to Assess Mating Potential in Space and Time
mathgraph	Directed and undirected graphs
matie	Measuring Association and Testing Independence Efficiently
matlab	MATLAB emulation package
matlabr	An Interface for MATLAB using System Calls
matlib	Matrix Functions for Teaching and Learning Linear Algebra and Multivariate Statistics
matpow	matrix powers
matR	Metagenomics Analysis Tools for R
Matrix	Sparse and Dense Matrix Classes and Methods
matrixcalc	Collection of functions for matrix calculations
MatrixCorrelation	Matrix Correlation Coefficients
MatrixEQTL	Matrix eQTL: Ultra fast eQTL analysis via large matrix operations
MatrixModels	Modelling with Sparse And Dense Matrices
matrixpls	Matrix-Based Partial Least Squares Estimation
matrixStats	Functions that Apply to Rows and Columns of Matrices (and to Vectors)
Matrix.utils	Data.frame-Like Operations on Sparse and Dense Matrix Objects
MATTOOLS	Modern Calibration Functions for the Modern Analog Technique (MAT)
MAVIS	Meta Analysis via Shiny
MAVTgsa	Three methods to identify differentially expressed gene sets, ordinary least square test, Multivariate Analysis Of Variance test with n contrasts and Random forest
MaXact	Exact max-type Cochran-Armitage trend test(CATT)
maxent	Low-memory Multinomial Logistic Regression with Support for Text Classification
MaxentVariableSelection	Selecting the Best Set of Relevant Environmental Variables along with the Optimal Regularization Multiplier for Maxent Niche Modeling
maxLik	Maximum Likelihood Estimation and Related Tools
maxlike	Model species distributions by estimating the probability of occurrence using presence-only data
MaxPro	Maximum Projection Designs
maxstat	Maximally Selected Rank Statistics
MazamaSpatialUtils	Mazama Science Spatial Data Download and Utility Functions
MBA	Multilevel B-spline Approximation
mbbefd	Maxwell Boltzmann Bose Einstein Fermi Dirac Distribution and Destruction Rate Modelling
MBCluster.Seq	Model-Based Clustering for RNA-seq Data
MBESS	The MBESS R Package
mbest	Moment-Based Estimation for Hierarchical Models
MBI	(M)ultiple-site (B)iodiversity (I)ndices Calculator
mblm	Median-Based Linear Models
MBmca	Nucleic Acid Melting Curve Analysis on Microbead Surfaces with R
mbmdr	Model Based Multifactor Dimensionality Reduction
mboost	Model-Based Boosting
MBTAr	Access Data from the Massachusetts Bay Transit Authority (MBTA) Web API
mBvs	Multivariate Bayesian Variable Selection Method Exploiting Dependence among Outcomes
mc2d	Tools for Two-Dimensional Monte-Carlo Simulations
MC2toPath	Translates information from netcdf files with MC2 output into inter-PVT transitions
MCAPS	MCAPS data and results
mcbiopi	Matrix Computation Based Identification Of Prime Implicants
mcc	Moment Corrected Correlation
mcclust	Process an MCMC Sample of Clusterings
MCDA	Functions to Support the Multicriteria Decision Aiding Process
MCDM	Multi-Criteria Decision Making Methods
mcemGLM	Maximum Likelihood Estimation for Generalized Linear Mixed Models
mcga	Machine Coded Genetic Algorithms for Real-Valued Optimization Problems
mcgibbsit	Warnes and Raftery's MCGibbsit MCMC diagnostic
mcGlobaloptim	Global optimization using Monte Carlo and Quasi Monte Carlo simulation
mcheatmaps	Multiple matrices heatmap visualization
MChtest	Monte Carlo hypothesis tests with Sequential Stopping
MCI	Multiplicative Competitive Interaction (MCI) Model
mcIRT	IRT models for multiple choice items (mcIRT)
MCL	Markov Cluster Algorithm
mcll	Monte Carlo Local Likelihood Estimation
mclogit	Mixed Conditional Logit
mclust	Normal Mixture Modelling for Model-Based Clustering, Classification, and Density Estimation
mcmc	Markov Chain Monte Carlo
MCMC4Extremes	Posterior Distribution of Extreme Value Models in R
MCMCglmm	MCMC Generalised Linear Mixed Models
MCMC.OTU	Bayesian Analysis of Multivariate Counts Data in DNA Metabarcoding and Ecology
MCMCpack	Markov Chain Monte Carlo (MCMC) Package
mcmcplots	Create Plots from MCMC Output
MCMC.qpcr	Bayesian Analysis of qRT-PCR Data
mcmcse	Monte Carlo Standard Errors for MCMC
mco	Multiple Criteria Optimization Algorithms and Related Functions
Mcomp	Data from the M-competitions
MConjoint	Conjoint Analysis through Averaging of Multiple Analyses
MCPAN	Multiple Comparisons Using Normal Approximation
mcparallelDo	A Simplified Interface for Running Commands on Parallel Processes
MCPerm	A Monte Carlo permutation method for multiple test correlation
MCPMod	Design and Analysis of Dose-Finding Studies (see also DoseFinding package)
mcprofile	Multiple Contrast Profiles
mcr	Method Comparison Regression
MCS	Model Confidence Set Procedure
mcsm	Functions for Monte Carlo Methods with R
McSpatial	Nonparametric spatial data analysis
MCTM	Markov Chains Transition Matrices
md	Selecting Bandwidth for Kernel Density Estimator with Minimum Distance Method
mda	Mixture and Flexible Discriminant Analysis
mdatools	Multivariate Data Analysis for Chemometrics
mded	Measuring the Difference Between Two Empirical Distributions
mdhglm	Multivariate Double Hierarchical Generalized Linear Models
MDimNormn	Multi-Dimensional MA Normalization for Plate Effect
MDM	Multinomial Diversity Model
MDMR	Multivariate Distance Matrix Regression
MDplot	Visualizing Molecular Dynamics Analyses
MDPtoolbox	Markov Decision Processes toolbox
MDR	Detect gene-gene interactions using multifactor dimensionality reduction
mdscore	Improved Score Tests for Generalized Linear Models
mdsdt	Functions for Analysis of Data with General Recognition Theory
MDSGUI	A GUI for interactive MDS in R
MeanShift	Clustering via the Mean Shift Algorithm
measuRing	Detection and Control of Tree-Ring Widths on Scanned Image Sections
meboot	Maximum Entropy Bootstrap for Time Series
mederrRank	Bayesian Methods for Identifying the Most Harmful Medication Errors
medflex	Flexible Mediation Analysis Using Natural Effect Models
MediaK	Calculate MeDiA_K Distance
Mediana	Clinical Trial Simulations
mediation	Causal Mediation Analysis
medicalrisk	Medical Risk and Comorbidity Tools for ICD-9-CM Data
MedOr	Median Ordering Statistical R package
medSTC	A max-margin supervised Sparse Topical Coding Model
MEET	MEET: Motif Elements Estimation Toolkit
mefa	Multivariate Data Handling in Ecology and Biogeography
mefa4	Multivariate Data Handling with S4 Classes and Sparse Matrices
megaptera	MEGAPhylogeny Techniques in R
meifly	Interactive model exploration using GGobi
mem	Moving Epidemic Method R Package
memgene	Spatial pattern detection in genetic distance data using Moran's Eigenvector Maps
memisc	Tools for Management of Survey Data and the Presentation of Analysis Results
memoise	Memoisation of Functions
MEMSS	Data sets from Mixed-effects Models in S
memuse	Memory Estimation Utilities
MenuCollection	Collection of Configurable GTK+ Menus
MergeGUI	A GUI for Merging Datasets in R
merror	Accuracy and Precision of Measurements
merTools	Tools for Analyzing Mixed Effect Regression Models
MESS	Miscellaneous Esoteric Statistical Scripts
meta	General Package for Meta-Analysis
meta4diag	Meta-Analysis for Diagnostic Test Studies
MetABEL	Meta-analysis of genome-wide SNP association results
MetabolAnalyze	Probabilistic latent variable models for metabolomic data
metabolomics	Analysis of Metabolomics Data
metacom	Analysis of the "Elements of Metacommunity Structure"
metacor	Meta-analysis of correlation coefficients
MetaCycle	Evaluate Periodicity in Large Scale Data
MetaDE	MetaDE: Microarray meta-analysis for differentially expressed gene detection
metafolio	Metapopulation simulations for conserving salmon through portfolio optimization
metafor	Meta-Analysis Package for R
metafuse	Fused Lasso Approach in Regression Coefficient Clustering
metagear	Comprehensive Research Synthesis Tools for Systematic Reviews and Meta-Analysis
metagen	Inference in Meta Analysis and Meta Regression
metaheur	Metaheuristic Optimization Framework for Preprocessing Combinations
MetaLandSim	Landscape and Range Expansion Simulation
metaLik	Likelihood Inference in Meta-Analysis and Meta-Regression Models
metaMA	Meta-analysis for MicroArrays
metamisc	Diagnostic and prognostic meta analysis (metamisc)
metaMix	Bayesian Mixture Analysis for Metagenomic Community Profiling
metansue	Meta-Analysis of Studies with Non Statistically-Significant Unreported Effects
metap	Meta-Analysis of Significance Values
MetaPath	Perform the Meta-Analysis for Pathway Enrichment Analysis (MAPE)
MetaPCA	MetaPCA: Meta-analysis in the Dimension Reduction of Genomic data
metaplus	Robust Meta-Analysis and Meta-Regression
MetaQC	MetaQC: Objective Quality Control and Inclusion/Exclusion Criteria for Genomic Meta-Analysis
metaRNASeq	Meta-analysis of RNA-seq data
metaSEM	Meta-Analysis using Structural Equation Modeling
metasens	Advanced Statistical Methods to Model and Adjust for Bias in Meta-Analysis
MetaSKAT	Meta Analysis for SNP-Set (Sequence) Kernel Association Test
metatest	Fit and test metaregression models
Metatron	Meta-analysis for Classification Data and Correction to Imperfect Reference
meteo	Spatio-Temporal Analysis and Mapping of Meteorological Observations
meteoForecast	Numerical Weather Predictions
meteogRam	Tools for plotting meteograms
meteR	Fitting and Plotting Tools for the Maximum Entropy Theory of Ecology (METE)
MetFns	Analysis of Visual Meteor Data
Meth27QC	Meth27QC: sample quality analysis, and sample control analysis
MethComp	Functions for Analysis of Agreement in Method Comparison Studies
Methplot	Visualize the methylation patterns
MethylCapSig	Detection of Differentially Methylated Regions using MethylCap-Seq Data
MetNorm	Statistical Methods for Normalizing Metabolomics Data
Metrics	Evaluation metrics for machine learning
metricsgraphics	Create Interactive Charts with the JavaScript 'MetricsGraphics' Library
metRology	Support for metrological applications
mets	Analysis of Multivariate Event Times
MetSizeR	GUI Tool for Estimating Sample Sizes for Metabolomic Experiments
MetStaT	Statistical metabolomics tools
mev	Multivariate Extreme Value Distributions
mewAvg	A Fixed Memeory Moving Expanding Window Average
MExPosition	Multi-table ExPosition
MF	Mitigated Fraction
MFAg	Multiple Factor Analysis (MFA)
MFHD	Multivariate Functional Halfspace Depth
mFilter	Miscellaneous time series filters
mfp	Multivariable Fractional Polynomials
MfUSampler	Multivariate-from-Univariate (MfU) MCMC Sampler
mfx	Marginal Effects, Odds Ratios and Incidence Rate Ratios for GLMs
mgarchBEKK	Simulating, Estimating and Diagnosing MGARCH (BEKK and mGJR) Processes
mgcv	Mixed GAM Computation Vehicle with GCV/AIC/REML Smoothness Estimation
MGGM	Structural Pursuit Over Multiple Undirected Graphs
MGLM	Multivariate Response Generalized Linear Models
mglmn	Model Averaging for Multivariate GLM with Null Models
mgm	Estimating Mixed Graphical Models
mgpd	mgpd: Functions for multivariate generalized Pareto distribution (MGPD of Type II)
mgraph	Graphing map attributes and non-map variables in R
MGRASTer	API Client for the MG-RAST Server of the US DOE KBase
MGSDA	Multi-Group Sparse Discriminant Analysis
mGSZ	Gene set analysis based on GSZ-scoring function and asymptotic p-value
MHadaptive	General Markov Chain Monte Carlo for Bayesian Inference using adaptive Metropolis-Hastings sampling
mhde	Minimum Hellinger Distance Test for Normality
mHG	Minimum-Hypergeometric Test
mhsmm	Inference for Hidden Markov and Semi-Markov Models
mht	Multiple Hypothesis Testing for Variable Selection in High-Dimensional Linear Models
MHTrajectoryR	Bayesian Model Selection in Logistic Regression for the Detection of Adverse Drug Reactions
mhurdle	Multiple hurdle Tobit models
mi	Missing Data Imputation and Model Checking
mice	Multivariate Imputation by Chained Equations
miceadds	Some Additional Multiple Imputation Functions, Especially for 'mice'
micEcon	Microeconomic Analysis and Modelling
micEconAids	Demand Analysis with the Almost Ideal Demand System (AIDS)
micEconCES	Analysis with the Constant Elasticity of Substitution (CES) function
micEconSNQP	Symmetric Normalized Quadratic Profit Function
miCoPTCM	Promotion Time Cure Model with Mis-Measured Covariates
microbats	An Implementation of Bat Algorithm in R
microbenchmark	Accurate Timing Functions
MicroDatosEs	Utilities for Official Spanish Microdata
micromap	Linked Micromap Plots
micromapST	Linked Micromap Plots for U. S. States
micropan	Microbial Pan-genome Analysis
MicroStrategyR	MicroStrategyR Package
MicSim	Performing Continuous-Time Microsimulation
midasr	Mixed Data Sampling Regression
midastouch	Multiple Imputation by Distance Aided Donor Selection
midrangeMCP	Multiple Comparisons Procedures Based on Studentized Midrange and Range Distributions
MigClim	Implementing dispersal into species distribution models
migest	Methods for the Indirect Estimation of Bilateral Migration
migration.indices	Migration indices
migui	Graphical User Interface to the 'mi' Package
MIICD	Multiple Imputation for Interval Censored Data
MIIVsem	Model Implied Instrumental Variable (MIIV) Estimation of Structural Equation Models
MILC	MIcrosimulation Lung Cancer (MILC) model
mime	Map Filenames to MIME Types
MImix	Mixture summary method for multiple imputation
MindOnStats	Data sets included in Utts and Heckard's Mind on Statistics
minerva	Maximal Information-Based Nonparametric Exploration R Package for Variable Analysis
Miney	Implementation of the Well-Known Game to Clear Bombs from a Given Field (Matrix)
miniCRAN	Create a Mini Version of CRAN Containing Only Selected Packages
miniGUI	tktcl quick and simple function GUI
minimap	Create Tile Grid Maps
minimax	Minimax distribution family
minimist	Parse Argument Options
miniUI	Shiny UI Widgets for Small Screens
minpack.lm	R Interface to the Levenberg-Marquardt Nonlinear Least-Squares Algorithm Found in MINPACK, Plus Support for Bounds
minPtest	Gene region-level testing procedure for SNP data, using the min P test resampling approach
minqa	Derivative-free optimization algorithms by quadratic approximation
minque	An R Package for Linear Mixed Model Analyses
MInt	Learn Direct Interaction Networks
minval	MINimal VALidation for Stoichiometric Reactions
minxent	Entropy Optimization Distributions
mipfp	Multidimensional Iterative Proportional Fitting and Alternative Models
MIPHENO	Mutant Identification through Probabilistic High throughput Enabled NOrmalization
miRada	MicroRNA Microarray Data Analysis
MiRSEA	'MicroRNA' Set Enrichment Analysis
mirt	Multidimensional Item Response Theory
mirtCAT	Computerized Adaptive Testing with Multidimensional Item Response Theory
miRtest	combined miRNA- and mRNA-testing
misc3d	Miscellaneous 3D Plots
miscF	Miscellaneous Functions
miscFuncs	Miscellaneous Useful Functions
miscset	Miscellaneous Tools Set
miscTools	Miscellaneous Tools and Utilities
MiSPU	Microbiome Based Sum of Powered Score (MiSPU) Tests
missDeaths	'Correctly Analyse Disease Recurrence with Missing at Risk Information using Population Mortality'
missForest	Nonparametric Missing Value Imputation using Random Forest
MissingDataGUI	A GUI for Missing Data Exploration
missMDA	Handling Missing Values with Multivariate Data Analysis
MissMech	Testing Homoscedasticity, Multivariate Normality, and Missing Completely at Random
MiST	Mixed effects Score Test for continuous outcomes
mistat	Data Sets, Functions and Examples from the Book: "Modern Industrial Statistics" by Kenett, Zacks and Amberti
mistral	Methods in Structural Reliability
MitISEM	Mixture of Student t Distributions using Importance Sampling and Expectation Maximization
mitml	Tools for Multiple Imputation in Multilevel Modeling
mitools	Tools for multiple imputation of missing data
mix	Estimation/Multiple Imputation for Mixed Categorical and Continuous Data
mixAK	Multivariate Normal Mixture Models and Mixtures of Generalized Linear Mixed Models Including Model Based Clustering
MixAll	Clustering using Mixture Models
mixcat	Mixed effects cumulative link and logistic regression models
mixdist	Finite Mixture Distribution Models
MixedDataImpute	Missing Data Imputation for Continuous and Categorical Data using Nonparametric Bayesian Joint Models
mixedMem	Tools for Discrete Multivariate Mixed Membership Models
MixedPoisson	Mixed Poisson Models
MixedTS	Mixed Tempered Stable Distribution
mixer	Random graph clustering
mixexp	Design and Analysis of Mixture Experiments
MIXFIM	Evaluation of the FIM in NLMEMs using MCMC
MixGHD	Model Based Clustering, Classification and Discriminant Analysis Using the Mixture of Generalized Hyperbolic Distributions
mixlm	Mixed Model ANOVA and Statistics for Education
MixMAP	Implements the MixMAP Algorithm
mixOmics	Omics Data Integration Project
mixor	Mixed-Effects Ordinal Regression Analysis
mixpack	Tools to Work with Mixture Components
mixPHM	Mixtures of Proportional Hazard Models
mixRasch	Mixture Rasch Models with JMLE
mixreg	Functions to fit mixtures of regressions
mixsep	Forensic Genetics DNA Mixture Separation
MixSim	Simulating Data to Study Performance of Clustering Algorithms
mixsmsn	Fitting Finite Mixture of Scale Mixture of Skew-Normal Distributions
mixtNB	DE Analysis of RNA-Seq Data by Mixtures of NB
mixtools	Tools for Analyzing Finite Mixture Models
mixtox	Curve Fitting and Mixture Toxicity Assessment
mixture	Mixture Models for Clustering and Classification
MixtureInf	Inference for Finite Mixture Models
mizer	Multi-species sIZE spectrum modelling in R
mkde	2D and 3D movement-based kernel density estimates (MKDEs)
mkin	Routines for Fitting Kinetic Models with One or More State Variables to Chemical Degradation Data
MKLE	Maximum kernel likelihood estimation
MKmisc	Miscellaneous Functions from M. Kohl
mkssd	Efficient multi-level k-circulant supersaturated designs
mlbench	Machine Learning Benchmark Problems
MLCIRTwithin	Latent Class Item Response Theory (LC-IRT) Models under Within-Item Multidimensionality
MLCM	Maximum Likelihood Conjoint Measurement
mlDNA	Machine Learning-based Differential Network Analysis of Transcriptome Data
mldr	Exploratory Data Analysis and Manipulation of Multi-Label Data Sets
mldr.datasets	R Ultimate Multilabel Dataset Repository
MLDS	Maximum Likelihood Difference Scaling
mlearning	Machine learning algorithms with unified interface and confusion matrices
MLEcens	Computation of the MLE for bivariate (interval) censored data
mlegp	Maximum Likelihood Estimates of Gaussian Processes
mleur	Maximum likelihood unit root test
mlgt	Multi-Locus Geno-Typing
mlica2	Independent Component Analysis using Maximum Likelihood
mlma	Multilevel Mediation Analysis
MLmetrics	Machine Learning Evaluation Metrics
mlmmm	ML estimation under multivariate linear mixed models with missing values
mlmRev	Examples from Multilevel Modelling Software Review
mlogit	multinomial logit model
mlogitBMA	Bayesian Model Averaging for Multinomial Logit Models
mlPhaser	Multi-Locus Haplotype Phasing
mlr	Machine Learning in R
MLRMPA	A package for Multilinear Regression Model Population Analysis
mlsjunkgen	Use the MLS Junk Generator Algorithm to Generate a Stream of Pseudo-Random Numbers
mlt	Most Likely Transformations
mlt.docreg	Most Likely Transformations: Documentation and Regression Tests
mlVAR	Multi-Level Vector Autoregression
mlxR	Simulation of Longitudinal Data
MM	The multiplicative multinomial distribution
MM2S	Single-Sample Classifier of Medulloblastoma Subtypes for Medulloblastoma Patient Samples, Mouse Models, and Cell Lines
MM2Sdata	Gene Expression Datasets for the 'MM2S' Package
mma	Multiple Mediation Analysis
mmand	Mathematical Morphology in Any Number of Dimensions
mmap	Map Pages of Memory
mmc	Multivariate Measurement Error Correction
mmcm	Modified Maximum Contrast Method
mmds	Mixture Model Distance Sampling (mmds)
mme	Multinomial Mixed Effects Models
mmeln	Estimation of Multinormal Mixture Distribution
mmeta	Multivariate Meta-Analysis
mmm	an R package for analyzing multivariate longitudinal data with multivariate marginal models
mmm2	Multivariate marginal models with shared regression parameters
MMMS	Multi-Marker Molecular Signature for Treatment-specific Subgroup Identification
mmod	Modern Measures of Population Differentiation
mmpp	Various Similarity and Distance Metrics for Marked Point Processes
mmppr	Markov Modulated Poisson Process for Unsupervised Event Detection in Time Series of Counts
MMS	Fixed effects Selection in Linear Mixed Models
mmtfa	Model-Based Clustering and Classification with Mixtures of Modified t Factor Analyzers
MMWRweek	Convert Dates to MMWR Day, Week, and Year
mnlogit	Multinomial Logit Model
MNM	Multivariate Nonparametric Methods. An Approach Based on Spatial Signs and Ranks
mnormpow	Multivariate Normal Distributions with Power Integrand
mnormt	The Multivariate Normal and t Distributions
MNP	R Package for Fitting the Multinomial Probit Model
MNS	Mixed Neighbourhood Selection
Mobilize	Mobilize plots and functions
MOCCA	Multi-objective optimization for collecting cluster alternatives
Modalclust	Hierarchical Modal Clustering
modeest	Mode Estimation
modehunt	Multiscale Analysis for Density Functions
modelfree	Model-free estimation of a psychometric function
ModelGood	Validation of risk prediction models
ModelMap	Modeling and Map Production using Random Forest and Stochastic Gradient Boosting
modelObj	A Model Object Framework for Regression Analysis
modeltools	Tools and Classes for Statistical Models
modes	Find the Modes and Assess the Modality of Complex and Mixture Distributions, Especially with Big Datasets
modiscloud	R tools for processing Level 2 Cloud Mask products from MODIS
MODISTools	MODIS Subsetting Tools
modMax	Community Structure Detection via Modularity Maximization
modQR	Multiple-Output Directional Quantile Regression
modTempEff	Modelling temperature effects using time series data
moduleColor	Basic Module Functions
modules	Self Contained Units of Source Code
mogavs	Multiobjective Genetic Algorithm for Variable Selection in Regression
MOJOV	Mojo Variants: Rare Variants analysis
mokken	Mokken Scale Analysis in R
molaR	Dental Surface Complexity Measurement Tools
mombf	Moment and Inverse Moment Bayes Factors
momentchi2	Moment-Matching Methods for Weighted Sums of Chi-Squared Random Variables
moments	Moments, cumulants, skewness, kurtosis and related tests
Momocs	Morphometrics using R
momr	Mining Metaomics Data (MetaOMineR)
mondate	Keep track of dates in terms of months
Mondrian	A Simple Graphical Representation of the Relative Occurrence and Co-Occurrence of Events
MonetDBLite	In-Process Version of MonetDB for R
MonetDB.R	Connect MonetDB to R
mongolite	Fast and Simple MongoDB Client for R
monitoR	Acoustic Template Detection in R
monmlp	Monotone Multi-Layer Perceptron Neural Network
monogeneaGM	Geometric Morphometric Analysis of Monogenean Anchors
monographaR	Taxonomic Monographs Tools
monomvn	Estimation for Multivariate Normal and Student-t Data with Monotone Missingness
MonoPhy	Allows to Explore Monophyly (or Lack of it) of Taxonomic Groups in a Phylogeny
MonoPoly	Functions to Fit Monotone Polynomials
monreg	Nonparametric Monotone Regression
moonBook	Functions and Datasets for the Book by Keon-Woong Moon
moonsun	Basic astronomical calculations with R
mopsocd	MOPSOCD: Multi-objective Particle Swarm Optimization with Crowding Distance
MOrder	Check Time Homogeneity and Markov Chain Order
morgenstemning	Color schemes compatible with red-green color perception difficulties
Morpho	Calculations and Visualisations Related to Geometric Morphometrics
morse	MOdelling Tools for Reproduction and Survival Data in Ecotoxicology
MorseGen	Simple raw data generator based on user-specified summary statistics
MortalitySmooth	Smoothing and Forecasting Poisson Counts with P-Splines
mosaic	Project MOSAIC Statistics and Mathematics Teaching Utilities
mosaicData	Project MOSAIC (mosaic-web.org) data sets
MoTBFs	Learning Hybrid Bayesian Networks using Mixtures of Truncated Basis Functions
MotilityLab	Quantitative Analysis of Motion
moult	Models for analysing moult in birds
mountainplot	Mountain Plots, Folded Empirical Cumulative Distribution Plots
mousetrack	Mouse-Tracking Measures from Trajectory Data
mousetrap	Process and Analyze Mouse-Tracking Data
move	Visualizing and Analyzing Animal Track Data
moveHMM	Animal Movement Modelling using Hidden Markov Models
movMF	Mixtures of von Mises-Fisher Distributions
mp	Multidimensional Projection Techniques
mpa	CoWords Method
MPAgenomics	Multi-Patient Analysis of Genomic Markers
mpath	Regularized Linear Models
mpbart	Multinomial Probit Bayesian Additive Regression Trees
MPCI	Multivariate Process Capability Indices (MPCI)
mpcv	Multivariate Process Capability Vector
MPDiR	Data sets and scripts for Modeling Psychophysical Data in R
mph	Multiscale persistent homology
MPINet	The package can implement the network-based metabolite pathway identification of pathways
MPLikelihoodWB	Modified Profile Likelihood Estimation for Weibull Shape and Regression Parameters
mplot	Graphical Model Stability and Variable Selection Procedures
MplusAutomation	Automating Mplus Model Estimation and Interpretation
mpm	Multivariate Projection Methods
mpMap	Multi-parent RIL genetic analysis
mpmcorrelogram	Multivariate Partial Mantel Correlogram
mpmi	Mixed-pair mutual information estimators
mpoly	Symbolic Computation and More with Multivariate Polynomials
Mposterior	Mposterior: R package for Robust and Scalable Bayes via a Median of Subset Posterior Measures
mppa	Statistics for analysing multiple simultaneous point processes on the real line
MPSEM	Modeling Phylogenetic Signals using Eigenvector Maps
mpt	Multinomial Processing Tree Models
MPTinR	Analyze Multinomial Processing Tree Models
mptools	RAMAS Metapop Tools
MPV	Data Sets from Montgomery, Peck and Vining's Book
mQTL	Metabolomic Quantitative Trait Locus Mapping
mra	Analysis of Mark-Recapture Data
mratios	Inferences for ratios of coefficients in the general linear model
MRCE	Multivariate regression with covariance estimation
MRCV	Methods for Analyzing Multiple Response Categorical Variables (MRCVs)
mrds	Mark-Recapture Distance Sampling
mreg	Fits regression models when the outcome is partially missing
MRH	Multi-Resolution Estimation of the Hazard Rate
mri	Modified Rand Index (1 and 2.1 and 2.2) and Modified Adjusted Rand Index (1 and 2.1 and 2.2)
MRIaggr	Management, Display, and Processing of Medical Imaging Data
mritc	MRI Tissue Classification
mRm	An R package for conditional maximum likelihood estimation in mixed Rasch models
mrMLM	Multilocus Random Mixed Linear Model for Genome-Wide Association Study
MRMR	Multivariate Regression Models for Reserving
mRMRe	R package for parallelized mRMR ensemble feature selection
MRQoL	Minimal Clinically Important Difference and Response Shift Effect for Health-Related Quality of Life
MRS	Multi-Resolution Scanning for Cross-Sample Differences
MRSP	Multinomial Response Models with Structured Penalties
MRsurv	A multiplicative-regression model for relative survival
MRwarping	Multiresolution time warping for functional data
msap	Statistical analysis for Methylation-sensitive Amplification Polymorphism data
msarc	Draw Diagrams (mis)Representing the Results of Mass Spec Experiments
msBP	Multiscale Bernstein Polynomials for Densities
MSBVAR	Markov-Switching, Bayesian, Vector Autoregression Models
MScombine	Combine Data from Positive and Negative Ionization Mode Finding Common Entities
msda	Multi-Class Sparse Discriminant Analysis
mseapca	Metabolite set enrichment analysis for factor loading in principal component analysis
MSeasy	Preprocessing of Gas Chromatography-Mass Spectrometry (GC-MS) data
MSeasyTkGUI	MSeasy Tcl/Tk Graphical User Interface
MSG	Data and Functions for the Book Modern Statistical Graphics
msgl	High Dimensional Multiclass Classification Using Sparse Group Lasso
msgpackR	A library to serialize or unserialize data in MessagePack format
msgps	Degrees of freedom of elastic net, adaptive lasso and generalized elastic net
msir	Model-based Sliced Inverse Regression
MSIseq	Assess Tumor Microsatellite Instability with a Decision Tree Classifier from Exome Somatic Mutations
msltrend	Improved Techniques to Estimate Trend, Velocity and Acceleration from Sea Level Records
msm	Multi-State Markov and Hidden Markov Models in Continuous Time
msma	Multiblock Sparse Multivariable Analysis
msme	Functions and Datasets for "Methods of Statistical Model Estimation"
msos	Datasets and Functions used in Multivariate Statistics: Old School by John Marden
MSQC	Multivariate Statistical Quality Control
msr	Morse-Smale Approximation, Regression and Visualization
ms.sev	Package for Calculation of ARMSS, Local MSSS and Global MSSS
msSurv	Nonparametric Estimation for Multistate Models
MST	Multivariate Survival Trees
mstate	Data Preparation, Estimation and Prediction in Multi-State Models
MSwM	Fitting Markov Switching Models
mtconnectR	Read Data from Delimited MTConnect Data Files
mtk	Mexico ToolKit library (MTK)
MTS	All-Purpose Toolkit for Analyzing Multivariate Time Series (MTS) and Estimating Multivariate Volatility Models
mtsdi	Multivariate time series data imputation
MTurkR	R Client for the MTurk Requester API
MTurkRGUI	A Graphical User Interface for MTurkR
MUCflights	Munich Franz-Josef-Strauss Airport Pattern Analysis
MuFiCokriging	Multi-Fidelity Cokriging models
muhaz	Hazard Function Estimation in Survival Analysis
muir	Exploring Data with Tree Data Structures
MultAlloc	Optimal Allocation in Stratified Sampling
multcomp	Simultaneous Inference in General Parametric Models
multcompView	Visualizations of Paired Comparisons
MultEq	Multiple Equivalence Tests and Simultaneous Confidence Intervals
multgee	GEE Solver for Correlated Nominal or Ordinal Multinomial Responses
multiAssetOptions	Finite Difference Method for Multi-Asset Option Valuation
multiband	Period Estimation for Multiple Bands
multibiplotGUI	Multibiplot Analysis in R
multic	Quantitative Linkage Analysis Tools using the Variance Components Approach
MultiCNVDetect	Multiple Copy Number Variation Detection
multicon	Multivariate Constructs
multicool	Permutations of Multisets in Cool-Lex Order
multifwf	Read Fixed Width Format Files Containing Lines of Different Type
MultiGHQuad	Multidimensional Gauss-Hermite Quadrature
multigroup	Multigroup Data Analysis
MultiLCIRT	Multidimensional Latent Class Item Response Theory Models
multilevel	Multilevel Functions
multilevelPSA	Multilevel Propensity Score Analysis
multimark	Capture-Mark-Recapture Analysis using Multiple Non-Invasive Marks
MultiMeta	Meta-analysis of Multivariate Genome Wide Association Studies
multinbmod	Regression analysis of overdispersed correlated count data
MultinomialCI	Simultaneous confidence intervals for multinomial proportions according to the method by Sison and Glaz
multinomRob	Robust Estimation of Overdispersed Multinomial Regression Models
MultiOrd	Generation of multivariate ordinal variates
MultiPhen	A Package to Test for Pleiotropic Effects
multiPIM	Variable Importance Analysis with Population Intervention Models
multipleNCC	Weighted Cox-regression for Nested Case-Control Data
multiplex	Algebraic Tools for the Analysis of Multiple Social Networks
multipol	multivariate polynomials
multirich	Calculate Multivariate Richness via UTC and sUTC
MultiRR	Bias, Precision, and Power for Multi-Level Random Regressions
multisensi	Multivariate Sensitivity Analysis
multisom	Clustering a Dataset using Multi-SOM Algorithm
multispatialCCM	Multispatial Convergent Cross Mapping
MultiSV	MultiSV: an R package for identification of structural variations in multiple populations based on whole genome resequencing
multitable	Simultaneous manipulation of multiple arrays of data, with data.list objects
multitaper	Multitaper Spectral Analysis Tools
MultivariateRandomForest	Multivariate Random Forest for Linearly Related Output Features
multivator	A multivariate emulator
multiwave	Estimation of Multivariate Long-Memory Models Parameters
multiway	Component Models for Multi-Way Data
multiwayvcov	Multi-way Standard Error Clustering
MultNonParam	Multivariate Nonparametric Methods
multxpert	Common Multiple Testing Procedures and Gatekeeping Procedures
muma	Metabolomics Univariate and Multivariate Analysis
MuMIn	Multi-Model Inference
munfold	Metric Unfolding
munsell	Utilities for Using Munsell Colours
munsellinterpol	Interpolate Munsell Renotation Data from Hue/Chroma to CIE/sRGB
muRL	Mailmerge using R, LaTeX, and the Web
murphydiagram	Murphy Diagrams for Forecast Comparisons
musicNMR	Conversion of Nuclear Magnetic Resonance spectrum in audio file
muStat	Prentice Rank Sum Test and McNemar Test
mutoss	Unified Multiple Testing Procedures
mutossGUI	A Graphical User Interface for the MuToss Project
MuViCP	MultiClass Visualizable Classification using Combination of Projections
MVA	An Introduction to Applied Multivariate Analysis with R
mvabund	Statistical Methods for Analysing Multivariate Abundance Data
MVar.pt	Analise multivariada (brazilian portuguese)
MVB	Mutivariate Bernoulli log-linear model
MvBinary	Modelling Multivariate Binary Data with Blocks of Specific One-Factor Distribution
mvbutils	Workspace organization, code and documentation editing, package prep and editing, etc
mvc	Multi-View Clustering
mvctm	Multivariate Variance Components Tests for Multilevel Data
mvcwt	Wavelet analysis of multiple time series
mvglmmRank	Multivariate Generalized Linear Mixed Models for Ranking Sports Teams
mvinfluence	Influence Measures and Diagnostic Plots for Multivariate Linear Models
mvmesh	Multivariate Meshes and Histograms in Arbitrary Dimensions
mvmeta	Multivariate and Univariate Meta-Analysis and Meta-Regression
mvMORPH	Multivariate Comparative Tools for Fitting Evolutionary Models to Morphometric Data
MVN	Multivariate Normality Tests
mvna	Nelson-Aalen estimator of the cumulative hazard in multistate models
mvnfast	Fast Multivariate Normal Methods
mvngGrAd	Moving Grid Adjustment in Plant Breeding Field Trials
mvnmle	ML estimation for multivariate normal data with missing values
mvnormtest	Normality test for multivariate variables
mvnpermute	Generate New Multivariate Normal Samples from Permutations
mvnTest	Goodness of Fit Tests for Multivariate Normality
mvoutlier	Multivariate outlier detection based on robust methods
mvProbit	Multivariate Probit Models
mvprpb	Orthant Probability of the Multivariate Normal Distribution
mvQuad	Methods for Multivariate Quadrature
MVR	Mean-Variance Regularization
mvrtn	Mean and Variance of Truncated Normal Distribution
mvsf	Shapiro-Francia Multivariate Normality Test
mvShapiroTest	Generalized Shapiro-Wilk test for multivariate normality
mvSLOUCH	Multivariate Stochastic Linear Ornstein-Uhlenbeck Models for Phylogenetic Comparative Hypotheses
MVT	Estimation and Testing for the Multivariate t-Distribution
mvtboost	Tree Boosting for Multivariate Outcomes
mvtmeta	Multivariate meta-analysis
mvtnorm	Multivariate Normal and t Distributions
mvtsplot	Multivariate Time Series Plot
mwa	Causal Inference in Spatiotemporal Event Data
mwaved	Multichannel Wavelet Deconvolution with Additive Long Memory Noise
mxkssd	Efficient mixed-level k-circulant supersaturated designs
MXM	Discovering Multiple, Statistically-Equivalent Signatures
mycobacrvR	Integrative immunoinformatics for Mycobacterial diseases in R platform
mycor	Automatic Correlation and Regression Test in a Data Frame
myepisodes	MyEpisodes RSS/API functions
Myrrix	Interface to Myrrix. Myrrix is a complete, real-time, scalable clustering and recommender system, evolved from Apache Mahout
Myrrixjars	External jars required for package Myrrix
myTAI	Performing Phylotranscriptomics with R
mztwinreg	Regression Models for Monozygotic Twin Data
0	
nabor	Wraps 'libnabo', a Fast K Nearest Neighbour Library for Low Dimensions
NADA	Nondetects And Data Analysis for environmental data
nadiv	(Non)Additive Genetic Relatedness Matrices
NAM	Nested Association Mapping Analysis
namespace	Provide namespace managment functions not (yet) present in base R
nanop	Tools for Nanoparticle Simulation and Calculation of PDF and Total Scattering Structure Function
NanoStringNorm	Normalize NanoString miRNA and mRNA Data
NAPPA	Performs the Processing and Normalisation of Nanostring miRNA and mRNA Data
nasaweather	Collection of datasets from the ASA 2006 data expo
nat	NeuroAnatomy Toolbox for Analysis of 3D Image Data
nat.nblast	NeuroAnatomy Toolbox (nat) Extension for Assessing Neuron Similarity and Clustering
nat.templatebrains	NeuroAnatomy Toolbox ('nat') Extension for Handling Template Brains
naturalsort	Natural Ordering
nat.utils	File System Utility Functions for 'NeuroAnatomy Toolbox'
NB	Maximum Likelihood method in estimating effective population size from genetic data
NbClust	Determining the Best Number of Clusters in a Data Set
nbconvertR	Vignette Engine Wrapping IPython Notebooks
NBDdirichlet	NBD-Dirichlet Model of Consumer Buying Behavior for Marketing Research
nbpMatching	Functions for Optimal Non-Bipartite Matching
NBPSeq	Negative Binomial Models for RNA-Sequencing Data
NCA	Necessary Condition Analysis
nCal	Nonlinear Calibration
ncappc	NCA Calculation and Population PK Model Diagnosis
ncbit	retrieve and build NBCI taxonomic data
ncdf4	Interface to Unidata netCDF (Version 4 or Earlier) Format Data Files
ncdf4.helpers	Helper functions for use with the ncdf4 package
ncdf.tools	Easier 'NetCDF' File Handling
nCDunnett	Noncentral Dunnett's Test Distribution
ncf	Spatial Nonparametric Covariance Functions
ncg	Computes the noncentral gamma function
NCmisc	Miscellaneous Functions for Creating Adaptive Functions and Scripts
ncvreg	Regularization Paths for SCAD and MCP Penalized Regression Models
ndl	Naive Discriminative Learning
ndtv	Network Dynamic Temporal Visualizations
neariso	Near-Isotonic Regression
NeatMap	Non-clustered heatmap alternatives
needs	Attaches and Installs Packages
needy	needy
NEff	Calculating Effective Sizes Based on Known Demographic Parameters of a Population
negenes	Estimating the Number of Essential Genes in a Genome
neldermead	R port of the Scilab neldermead module
neotoma	Access to the Neotoma Paleoecological Database Through R
nephro	Biostatistics Utilities for Nephrology
NEpiC	Network Assisted Algorithm for Epigenetic Studies Using Mean and Variance Combined Signals
NestedCohort	Survival Analysis for Cohorts with Missing Covariate Information
nestedRanksTest	Mann-Whitney-Wilcoxon Test for Nested Ranks
netassoc	Inference of Species Associations from Co-Occurrence Data
netClass	netClass: An R Package for Network-Based Biomarker Discovery
NetCluster	Clustering for networks
netcoh	Statistical Modeling with Network Cohesion
NetComp	Network Generation and Comparison
NetData	Network Data for McFarland's SNA R labs
netdiffuseR	Network Analysis for Diffusion of Innovations
netgen	Network Generator for Combinatorial Graph Problems
netgsa	Network-Based Gene Set Analysis
NetIndices	Estimating network indices, including trophic structure of foodwebs in R
netmeta	Network Meta-Analysis using Frequentist Methods
NetPreProc	Network Pre-Processing and Normalization
nets	Network Estimation for Time Series
NetSim	A Social Networks Simulation Tool in R
NetSwan	Network Strengths and Weaknesses Analysis
nettools	A Network Comparison Framework
network	Classes for Relational Data
networkD3	D3 JavaScript Network Graphs from R
networkDynamic	Dynamic Extensions for Network Objects
networkDynamicData	Dynamic (Longitudinal) Network Datasets
networkreporting	Tools for using network reporting estimators
networksis	Simulate Bipartite Graphs with Fixed Marginals Through Sequential Importance Sampling
networkTomography	Tools for network tomography
neural	Neural Networks
neuralnet	Training of neural networks
NeuralNetTools	Visualization and Analysis Tools for Neural Networks
neuroblastoma	Neuroblastoma copy number profiles
neuroim	Data Structures and Handling for Neuroimaging Data
neuRosim	Functions to Generate fMRI Data Including Activated Data, Noise Data and Resting State Data
Newdistns	Computes Pdf, Cdf, Quantile and Random Numbers, Measures of Inference for 19 General Families of Distributions
nFactors	Parallel Analysis and Non Graphical Solutions to the Cattell Scree Test
nFCA	Numerical Formal Concept Analysis for Systematic Clustering
ngram	An n-gram Babbler
ngramrr	A Simple General Purpose N-Gram Tokenizer
ngspatial	Fitting the centered autologistic and sparse spatial generalized linear mixed models for areal data
NHANES	Data from the US National Health and Nutrition Examination Study
nhanesA	NHANES Data Retrieval
NHEMOtree	Non-hierarchical evolutionary multi-objective tree learner to perform cost-sensitive classification
nhlscrapr	Compiling the NHL Real Time Scoring System Database for easy use in R
NHMM	Bayesian NHMM Modeling (Multiple Time Series)
NHMSAR	Non-Homogeneous Markov Switching Autoregressive Models
NHPoisson	Modelling and Validation of Non Homogeneous Poisson Processes
nice	Get or Set UNIX Niceness
nicheROVER	(Niche) (R)egion and Niche (Over)lap Metrics for Multidimensional Ecological Niches
NightDay	Night and Day Boundary Plot Funtion
nima	Nima Hejazi's Miscellaneous R Code
nimble	Flexible BUGS-Compatible System for Hierarchical Statistical Modeling and Algorithm Development
Nippon	Japanese Utility Functions and Data
NIPTeR	Fast and Accurate Trisomy Prediction in Non-Invasive Prenatal Testing
NISTnls	Nonlinear least squares examples from NIST
NISTunits	Fundamental Physical Constants and Unit Conversions from NIST
nivm	Noninferiority Tests with Variable Margins
NlcOptim	Solve Nonlinear Optimization with Nonlinear Constraints
nleqslv	Solve Systems of Nonlinear Equations
nlme	Linear and Nonlinear Mixed Effects Models
nlmeODE	Non-linear mixed-effects modelling in nlme using differential equations
nlmeU	Datasets and utility functions enhancing functionality of nlme package
nlmrt	Functions for Nonlinear Least Squares Solutions
nlnet	Nonlinear Network Reconstruction and Clustering Based on DCOL (Distance Based on Conditional Ordered List)
nloptr	R interface to NLopt
NLP	Natural Language Processing Infrastructure
NLPutils	Natural Language Processing Utilities
nlreg	Higher Order Inference for Nonlinear Heteroscedastic Models
NLRoot	searching for the root of equation
nlrr	Non-Linear Relative Risk Estimation and Plotting
nls2	Non-linear regression with brute force
nlsem	Fitting Structural Equation Mixture Models
nlsMicrobio	Nonlinear regression in predictive microbiology
nlsmsn	Fitting nonlinear models with scale mixture of skew-normal distributions
nlstools	Tools for Nonlinear Regression Analysis
NlsyLinks	Utilities and Kinship Information for Research with the NLSY
nlt	A nondecimated lifting transform for signal denoising
nlts	(non)linear time series analysis
nLTT	Calculate the NLTT Statistic
nlWaldTest	Wald Test of Nonlinear Restrictions and Nonlinear CI
nmcdr	Non-parametric Multiple Change-points Detection
NMF	Algorithms and Framework for Nonnegative Matrix Factorization (NMF)
NMFN	Non-negative Matrix Factorization
NMOF	Numerical Methods and Optimization in Finance
nnet	Feed-Forward Neural Networks and Multinomial Log-Linear Models
nnetpredint	Prediction Intervals of Multi-Layer Neural Networks
nnlasso	Non-Negative Lasso and Elastic Net Penalized Generalized Linear Models
NNLM	Fast and Versatile Non-Negative Matrix Factorization
nnls	The Lawson-Hanson algorithm for non-negative least squares (NNLS)
NNTbiomarker	Calculate Design Parameters for Biomarker Validation Studies
nodeHarvest	Node Harvest for Regression and Classification
nodiv	Compares the Distribution of Sister Clades Through a Phylogeny
noia	Implementation of the Natural and Orthogonal InterAction (NOIA) model
nomclust	Hierarchical Nominal Clustering Package
NominalLogisticBiplot	Biplot representations of categorical data
noncensus	U.S. Census Regional and Demographic Data
noncompliance	Causal Inference in the Presence of Treatment Noncompliance Under the Binary Instrumental Variable Model
nonlinearTseries	Nonlinear Time Series Analysis
nonnest2	Tests of Non-Nested Models
nonparaeff	Nonparametric Methods for Measuring Efficiency and Productivity
NonpModelCheck	Model Checking and Variable Selection in Nonparametric Regression
nonrandom	Stratification and matching by the propensity score
nontarget	Detecting Isotope, Adduct and Homologue Relations in LC-MS Data
nontargetData	Quantized simulation data of isotope pattern centroids
nopp	Nash Optimal Party Positions
nor1mix	Normal (1-d) Mixture Models (S3 Classes and Methods)
nordklimdata1	Dataset for Climate Analysis with Data from the Nordic Region
norm	Analysis of multivariate normal datasets with missing values
NormalGamma	Normal-gamma convolution model
NormalLaplace	The Normal Laplace Distribution
normalp	Routines for Exponential Power Distribution
NormPsy	Normalisation of Psychometric Tests
NORMT3	Evaluates complex erf, erfc, Faddeeva, and density of sum of Gaussian and Student's t
normtest	Tests for Normality
normwhn.test	Normality and White Noise Testing
NORRRM	Geochemical Toolkit for R
NORTARA	Generation of Multivariate Data with Arbitrary Marginals
nortest	Tests for Normality
nose	nose Package for R
NostalgiR	Advanced Text-Based Plots
notifyR	Send push notifications to your smartphone via pushover.net (ACCOUNT REQUIRED!)
novelist	NOVEL Integration of the Sample and Thresholded (NOVELIST) Correlation and Covariance Estimators
noweb	Noweb system for R
Nozzle.R1	Nozzle Reports
np	Nonparametric kernel smoothing methods for mixed data types
nparACT	Non-Parametric Measures of Actigraphy Data
nparcomp	Multiple Comparisons and Simultaneous Confidence Intervals
nparLD	Nonparametric Analysis of Longitudinal Data in Factorial Experiments
NPBayesImpute	Non-Parametric Bayesian Multiple Imputation for Categorical Data
npbr	Nonparametric Boundary Regression
NPC	Nonparametric Combination of Hypothesis Tests
NPCD	Nonparametric Methods for Cognitive Diagnosis
NPCirc	Nonparametric Circular Methods
npcp	Some Nonparametric Tests for Change-Point Detection in Possibly Multivariate Observations
npde	Normalised prediction distribution errors for nonlinear mixed-effect models
NPHMC	Sample Size Calculation for the Proportional Hazards Mixture Cure Model
npIntFactRep	Nonparametric Interaction Tests for Factorial Designs with Repeated Measures
nplplot	Plotting linkage and association results
nplr	N-Parameter Logistic Regression
NPMLEcmprsk	Type-Specific Failure Rate and Hazard Rate on Competing Risks Data
npmlreg	Nonparametric maximum likelihood estimation for random effect models
NPMPM	tertiary probabilistic model in predictive microbiology for use in food manufacture
npmr	Nuclear Penalized Multinomial Regression
npmv	Nonparametric Comparison of Multivariate Samples
NPMVCP	Nonparametric Multivariate Change Point Model
nppbib	Nonparametric Partially-Balanced Incomplete Block Design Analysis
npregfast	Nonparametric Estimation of Regression Models with Factor-by-Curve Interactions
nproc	Neyman-Pearson Receiver Operator Curve
NPS	Convenience Functions and Tests for Working With the Net Promoter Score (NPS)
npsf	Nonparametric and Stochastic Efficiency and Productivity Analysis
NPsimex	Nonparametric Smoothing for contaminated data using Simulation-Extrapolation
npsm	Package for Nonparametric Statistical Methods using R
npsp	Nonparametric spatial (geo)statistics
npst	Generalization of Hewitt's Seasonality Test
npsurv	Non-Parametric Survival Analysis
NSA	Post-normalization of total copy numbers
nscancor	Non-Negative and Sparse CCA
NScluster	Simulation and Estimation of the Neyman-Scott Type Spatial Cluster Models
nsga2R	Elitist Non-dominated Sorting Genetic Algorithm based on R
nsgp	Non-Stationary Gaussian Process Regression
NSM3	Functions and Datasets to Accompany Hollander, Wolfe, and Chicken - Nonparametric Statistical Methods, Third Edition
nsprcomp	Non-Negative and Sparse PCA
nsRFA	Non-supervised Regional Frequency Analysis
NSUM	Network Scale Up Method
nullabor	Tools for Graphical Inference
numbers	Number-Theoretic Functions
numDeriv	Accurate Numerical Derivatives
numOSL	Numeric Routines for Optically Stimulated Luminescence Dating
nutshell	Data for "R in a Nutshell"
nutshell.audioscrobbler	Audioscrobbler data for "R in a Nutshell"
nutshell.bbdb	Baseball Database for "R in a Nutshell"
nws	R functions for NetWorkSpaces and Sleigh
nycflights13	Data about flights departing NYC in 2013
0	
oaColors	OpenAnalytics Colors Package
oai	General Purpose 'Oai-PMH' Services Client
OAIHarvester	Harvest Metadata Using OAI-PMH v2.0
oapackage	Orthogonal Array Package
oaPlots	OpenAnalytics Plots Package
Oarray	Arrays with arbitrary offsets
oasis	Multiple Sclerosis Lesion Segmentation using Magnetic Resonance Imaging (MRI)
OasisR	Outright Tool for the Analysis of Spatial Inequalities and Segregation
oaxaca	Blinder-Oaxaca Decomposition
obAnalytics	Limit Order Book Analytics
objectProperties	A factory of self-describing properties
objectSignals	objectSignals
obliclus	Cluster-based factor rotation
obliqueRF	Oblique Random Forests from Recursive Linear Model Splits
oblique.tree	Oblique Trees for Classification Data
obs.agree	An R package to assess agreement between observers
OBsMD	Objective Bayesian Model Discrimination in Follow-Up Designs
obsSens	Sensitivity analysis for Observational studies
oc	OC Roll Call Analysis Software
occ	Estimates PET neuroreceptor occupancies
oce	Analysis of Oceanographic Data
ocean	Biophysical Oceanography Tools
OceanView	Visualisation of Oceanographic Data and Model Output
ocedata	Oceanographic Datasets for Oce
ocomposition	Regression for Rank-Indexed Compositional Data
OData	R Helper for OData Web Services
ODB	Open Document Databases (.odb) management
odds.converter	Betting Odds Conversion
odeintr	C++ ODE Solvers Compiled on-Demand
odfWeave	Sweave processing of Open Document Format (ODF) files
odfWeave.survey	Support for odfWeave on the survey package
ODMconverter	Tools to Convert ODM Files
OECD	Search and Extract Data from the OECD
oem	Orthogonalizing Expectation maximization
oglmx	Estimation of Ordered Generalized Linear Models
Ohmage	R Client for Ohmage 2 server
OIdata	Data sets and supplements (OpenIntro)
OIsurv	Survival analysis supplement to OpenIntro guide
OjaNP	Multivariate Methods Based on the Oja Median and Related Concepts
okmesonet	Retrieve Oklahoma Mesonet climatological data
olctools	Open Location Code Handling in R
OligoSpecificitySystem	Oligo Specificity System
OLScurve	OLS growth curve trajectories
omd	filter the molecular descriptors for QSAR
OmicKriging	Poly-Omic Prediction of Complex TRaits
omics	'–omics' Data Analysis Toolbox
oncomodel	Maximum likelihood tree models for oncogenesis
Oncotree	Estimating oncogenetic trees
OneArmPhaseTwoStudy	Planing, Monitoring and Evaluating Oncological Phase 2 Studies
onemap	Software for constructing genetic maps in experimental crosses: full-sib, RILs, F2 and backcrosses
ONETr	Efficient Authenticated Interaction with the O*NET API
OneTwoSamples	Deal with one and two (normal) samples
onewaytests	One-Way Independent Groups Design
onion	octonions and quaternions
onlinePCA	Online Principal Component Analysis
onls	Orthogonal Nonlinear Least-Squares Regression
ontologyIndex	Functions for Reading Ontologies into R
ontologyPlot	Functions for Visualising Sets of Ontological Terms
ontologySimilarity	Functions for Calculating Ontological Similarities
OOmisc	Ozgur-Ozlem Miscellaneous
OpasnetUtils	Opasnet Modelling Environment Utility Functions
OPDOE	OPtimal Design Of Experiments
opefimor	Option Pricing and Estimation of Financial Models in R
openair	Tools for the Analysis of Air Pollution Data
OpenCL	Interface allowing R to use OpenCL
opencpu	A System for Embedded Scientific Computing and Reproducible Research with R
openintro	OpenIntro data sets and supplemental functions
OpenMPController	Control number of OpenMP threads dynamically
OpenMx	Extended Structural Equation Modelling
openNLP	Apache OpenNLP Tools Interface
openNLPdata	Apache OpenNLP Jars and Basic English Language Models
OpenRepGrid	Tools to analyse repertory grid data
openssl	Toolkit for Encryption, Signatures and Certificates Based on OpenSSL
OpenStreetMap	Access to Open Street Map Raster Images
opentraj	Tools for Creating and Analysing Air Trajectory Data
openVA	Automated Method for Verbal Autopsy
openxlsx	Read, Write and Edit XLSX Files
operators	Additional Binary Operators
operator.tools	Utilities for Working with R's Operators
OPI	Open Perimetry Interface
ops	Optimal Power Space Transformation
optAUC	Optimal Combinations of Diagnostic Tests Based on AUC
optBiomarker	Estimation of optimal number of biomarkers for two-group microarray based classifications at a given error tolerance level for various classification rules
optCluster	Determine Optimal Clustering Algorithm and Number of Clusters
optextras	Tools to Support Optimization Methods
OptGS	Near-Optimal and Balanced Group-Sequential Designs for Clinical Trials with Continuous Outcomes
OptHedging	Estimation of value and hedging strategy of call and put options
optifunset	Set Options if Unset
optigrab	Command-Line Parsing for an R World
OptimalCutpoints	Computing optimal cutpoints in diagnostic tests
OptimaRegion	Confidence Regions for Optima
optimbase	R port of the Scilab optimbase module
optimsimplex	R port of the Scilab optimsimplex module
optimx	A Replacement and Extension of the optim() Function
OptInterim	Optimal Two and Three Stage Designs for Single-Arm and Two-Arm Randomized Controlled Trials with a Long-Term Binary Endpoint
OptionPricing	Option Pricing with Efficient Simulation Algorithms
OptiQuantR	Simplifies and Automates Analyzing and Reporting OptiQuant's log Data
optiRum	Financial Functions & More
optiscale	Optimal scaling
optismixture	Optimal Mixture Weights in Multiple Importance Sampling
optmatch	Functions for Optimal Matching
optparse	Command Line Option Parser
optpart	Optimal Partitioning of Similarity Relations
optR	Optimization Toolbox for solving linear systems
optrees	Optimal Trees in Weighted Graphs
ora	Convenient Tools for Working with Oracle Databases
orca	Computation of Graphlet Orbit Counts in Sparse Graphs
ORCI	Several confidence intervals for the odds ratio
orclus	ORCLUS subspace clustering
ORCME	Order Restricted Clustering for Microarray Experiments
orcutt	Estimate procedure in case of first order autocorrelation
ordBTL	Modelling comparison data with ordinal response
orddom	Ordinal Dominance Statistics
ORDER2PARENT	Estimate parent distributions with data of several order statistics
orderbook	Orderbook visualization/Charting software
orderedLasso	Ordered Lasso and Time-lag Sparse Regression
OrdFacReg	Least Squares, Logistic, and Cox-Regression with Ordered Predictors
ordiBreadth	Ordinated Diet Breadth
ordinal	Regression Models for Ordinal Data
ordinalCont	Ordinal Regression Analysis for Continuous Scales
ordinalgmifs	Ordinal Regression for High-dimensional Data
OrdinalLogisticBiplot	Biplot representations of ordinal variables
ordinalNet	Penalized Ordinal Regression
OrdLogReg	Ordinal Logic Regression
OrdMonReg	Compute least squares estimates of one bounded or two ordered isotonic regression curves
OrdNor	Concurrent Generation of Ordinal and Normal Data with Given Correlation Matrix and Marginal Distributions
ordPens	Selection and/or Smoothing of Ordinal Predictors
ore	An R Interface to the Oniguruma Regular Expression Library
ores	Connector to the Objective Revision Evaluation Service (ORES)
OrgMassSpecR	Organic Mass Spectrometry
orgR	Analyse Text Files Created by Emacs' Org mode
ORIClust	Order-restricted Information Criterion-based Clustering Algorithm
orientlib	Support for orientation data
OriGen	Fast Spatial Ancestry via Flexible Allele Frequency Surfaces
orloca	The package deals with Operations Research LOCational Analysis models
orloca.es	Spanish version of orloca package
ORMDR	ORMDR
oro.dicom	Rigorous - DICOM Input / Output
oro.nifti	Rigorous - NIfTI + ANALYZE + AFNI : Input / Output
oro.pet	Rigorous - Positron Emission Tomography
orQA	Order Restricted Assessment Of Microarray Titration Experiments
orsifronts	Southern Ocean Frontal Distributions (Orsi)
orsk	Converting Odds Ratio to Relative Risk in Cohort Studies with Partial Data Information
orthogonalsplinebasis	Orthogonal B-Spline Basis Functions
OrthoPanels	Dynamic Panel Models with Orthogonal Reparameterization of Fixed Effects
orthopolynom	Collection of functions for orthogonal and orthonormal polynomials
osd	Orthogonal Signal Deconvolution for Spectra Deconvolution in GC-MS and GCxGC-MS Data
osDesign	Design and analysis of observational studies
osmar	OpenStreetMap and R
osmplotr	Customisable Images of OpenStreetMap Data
osrm	Interface Between R and the OpenStreetMap-Based Routing Service OSRM
OTE	Optimal Trees Ensembles for Regression, Classification and Class Membership Probability Estimation
OTRselect	Variable Selection for Optimal Treatment Decision
OTUtable	North Temperate Lakes - Microbial Observatory 16S Time Series Data and Functions
ouch	Ornstein-Uhlenbeck Models for Phylogenetic Comparative Hypotheses
outbreaker	Bayesian Reconstruction of Disease Outbreaks by Combining Epidemiologic and Genomic Data
OutbreakTools	Basic Tools for the Analysis of Disease Outbreaks
OutlierDC	Outlier Detection using quantile regression for Censored Data
OutlierDM	Outlier Detection for Multi-replicated High-throughput Data
outliers	Tests for outliers
OutrankingTools	Functions for Solving Multiple-criteria Decision-making Problems
OUwie	Analysis of Evolutionary Rates in an OU Framework
overlap	Estimates of Coefficient of Overlapping for Animal Activity Patterns
overlapping	Estimation of Overlapping in Empirical Distributions
OxyBS	Processing of Oxy-Bisulfite Microarray Data
oz	Plot the Australian Coastline and States
0	
P2C2M	Posterior Predictive Checks of Coalescent Models
p2distance	Welfare's Synthetic Indicator
p3state.msm	Analyzing survival data
pa	Performance Attribution for Equity Portfolios
PAactivPAL	Summarize Daily Physical Activity from 'activPAL' Accelerometer Data
PabonLasso	Pabon Lasso Graphs and Comparing Situations of a Unit in Two Different Times
PAC	Partition-Assisted Clustering
pacbpred	PAC-Bayesian Estimation and Prediction in Sparse Additive Models
pack	Convert values to/from raw vectors
packagetrackr	Track R Package Downloads from RStudio's CRAN Mirror
packcircles	Circle Packing
packClassic	Toy example of Pack Classic
packHV	A few useful functions for statisticians
packrat	A Dependency Management System for Projects and their R Package Dependencies
packS4	Toy Example of S4 Package
pacman	Package Management Tool
paco	Procrustes Application to Cophylogenetic Analysis
PACVB	Variational Bayes (VB) Approximation of Gibbs Posteriors with Hinge Losses
Pade	Padé Approximant Coefficients
paf	Attributable Fraction Function for Censored Survival Data
PAFit	Nonparametric Estimation of Preferential Attachment and Node Fitness in Temporal Complex Networks
pagenum	Put Page Numbers on Graphics
pageviews	An API Client for Wikimedia Traffic Data
PAGI	The package can identify the dysregulated KEGG pathways based on global influence from the internal effect of pathways and crosstalk between pathways
PAGWAS	Pathway Analysis Methods for Genomewide Association Data
pairedCI	Confidence intervals for the ratio of locations and for the ratio of scales of two paired samples
PairedData	Paired Data Analysis
pairheatmap	A tool for comparing heatmaps
pairsD3	D3 Scatterplot Matrices
PairViz	Visualization using Eulerian tours and Hamiltonian decompositions
pairwise	Rasch Model Parameters by Pairwise Algorithm
pairwiseCI	Confidence Intervals for Two Sample Comparisons
palaeoSig	Significance Tests for Palaeoenvironmental Reconstructions
paleobioDB	Download and Process Data from the Paleobiology Database
paleofire	Analysis of Charcoal Records from the Global Charcoal Database
paleoMAS	Paleoecological Analysis
paleotree	Paleontological and Phylogenetic Analyses of Evolution
paleoTS	Analyze Paleontological Time-Series
palettetown	Use Pokemon Inspired Colour Palettes
palinsol	Insolation for Palaeoclimate Studies
palr	Colour Palettes for Data
pamctdp	Principal Axes Methods for Contingency Tables with Partition Structures on Rows and Columns
pamm	Power Analysis for Random Effects in Mixed Models
pampe	Implementation of the Panel Data Approach Method for Program Evaluation
pamr	Pam: prediction analysis for microarrays
pan	Multiple Imputation for Multivariate Panel or Clustered Data
pAnalysis	Benchmarking and Rescaling R2 using Noise Percentile Analysis
PANDA	Preferential Attachment Based Common Neighbor Distribution Derived Functional Associations
pander	An R Pandoc Writer
panelaggregation	Aggregate Longitudinal Survey Data
panelAR	Estimation of Linear AR(1) Panel Data Models with Cross-Sectional Heteroskedasticity and/or Correlation
PanelCount	Random Effects and/or Sample Selection Models for Panel Count Data
Paneldata	Linear models for panel data
pangaear	Client for the 'Pangaea' Database
PANICr	PANIC Tests of Nonstationarity
papeR	A Toolbox for Writing Pretty Papers and Reports
ParallelForest	Random Forest Classification with Parallel Computing
parallelize.dynamic	Automate parallelization of function calls by means of dynamic code analysis
parallelMap	Unified Interface to Parallelization Back-Ends
parallelMCMCcombine	Methods for combining independent subset Markov chain Monte Carlo (MCMC) posterior samples to estimate a posterior density given the full data set
parallelML	A Parallel-Voting Algorithm for many Classifiers
ParallelPC	Paralellised Versions of Constraint Based Causal Discovery Algorithms
parallelSVM	A Parallel-Voting Version of the Support-Vector-Machine Algorithm
ParamHelpers	Helpers for Parameters in Black-Box Optimization, Tuning and Machine Learning
paramlink	Parametric Linkage Analysis in R
params	Simplify Parameters
paran	Horn's Test of Principal Components/Factors
parboost	Distributed Model-Based Boosting
parcor	Regularized estimation of partial correlation matrices
ParDNAcopy	Parallel implementation of the "segment" function of package "DNAcopy"
ParentOffspring	Conduct the Parent-Offspring Test Using Monomorphic SNP Markers
ParetoPosStable	Computing, Fitting and Validating the PPS Distribution
parfm	Parametric Frailty Models
parfossil	Parallelized functions for palaeoecological and palaeogeographical analysis
parma	Portfolio Allocation and Risk Management Applications
parmigene	Parallel Mutual Information estimation for Gene Network reconstruction
parsec	Partial Orders in Socio-Economics
parsedate	Recognize and Parse Dates in Various Formats, Including All ISO 8601 Formats
partDSA	Partitioning Using Deletion, Substitution, and Addition Moves
partialAR	Partial Autoregression
partialOR	Partial Odds Ratio
partitionMap	Partition Maps
partitionMetric	Compute a distance metric between two partitions of a set
partitions	Additive Partitions of Integers
partools	Tools for the 'Parallel' Package
partsm	Periodic Autoregressive Time Series Models
party	A Laboratory for Recursive Partytioning
partykit	A Toolkit for Recursive Partytioning
parviol	Parviol
PAS	Polygenic Analysis System (PAS)
Pasha	Preprocessing of Aligned Sequences from HTS Analyses
pass	Prediction and Stability Selection of Tuning Parameters
password	Create Random Passwords
pastecs	Package for Analysis of Space-Time Ecological Series
pastis	Phylogenetic Assembly with Soft Taxonomic Inferences
PASWR	PROBABILITY and STATISTICS WITH R
PASWR2	Probability and Statistics with R, Second Edition
patchDVI	Package to Patch .dvi or .synctex Files
patchPlot	Scatterplots of image patches
patchSynctex	Communication Between Editor and Viewer for Literate Programs
pathClass	Classification using biological pathways as prior knowledge
pathdiagram	Basic functions for drawing path diagrams
pathmox	Pathmox Approach of Segmentation Trees in Partial Least Squares Path Modeling
pathological	Path Manipulation Utilities
patPRO	Visualizing Temporal Microbiome Data
PatternClass	Class-focused pattern metric comparisons using simulation
pauwels2014	Bayesian Experimental Design for Systems Biology
pavo	Perceptual Analysis, Visualization and Organization of Spectral Color Data in R
pawacc	Physical activity with accelerometers
PAWL	Implementation of the PAWL algorithm
pbapply	Adding Progress Bar to '*apply' Functions
pbatR	P2BAT
PBD	Protracted Birth-Death Model of Diversification
pbdBASE	Programming with Big Data – Base Wrappers for Distributed Matrices
pbdDEMO	Programming with Big Data – Demonstrations and Examples Using 'pbdR' Packages
pbdDMAT	Programming with Big Data – Distributed Matrix Methods
pbdMPI	Programming with Big Data – Interface to MPI
pbdNCDF4	Programming with Big Data – Interface to Parallel Unidata NetCDF4 Format Data Files
pbdPROF	Programming with Big Data — MPI Profiling Tools
pbdSLAP	Programming with Big Data – Scalable Linear Algebra Packages
pbdZMQ	Programming with Big Data – Interface to ZeroMQ
PBImisc	A Set of Datasets Used in My Classes or in the Book 'Modele Liniowe i Mieszane w R, Wraz z Przykladami w Analizie Danych'
pbivnorm	Vectorized Bivariate Normal CDF
pbkrtest	Parametric Bootstrap and Kenward Roger Based Methods for Mixed Model Comparison
pbo	Probability of Backtest Overfitting
pBrackets	Plot Brackets
pbs	Periodic B Splines
PBSadmb	ADMB for R Using Scripts or GUI
PBSddesolve	Solver for Delay Differential Equations
PBSmapping	Mapping Fisheries Data and Spatial Analysis Tools
PBSmodelling	GUI Tools Made Easy: Interact with Models and Explore Data
pca3d	Three Dimensional PCA Plots
PCA4TS	Segmenting Multiple Time Series by Contemporaneous Linear Transformation
pcaBootPlot	Create 2D Principal Component Plots with Bootstrapping
pcadapt	Principal Component Analysis for Outlier Detection
pcaL1	Three L1-Norm PCA Methods
pcalg	Methods for Graphical Models and Causal Inference
PCAmixdata	Multivariate Analysis of Mixed Data
pcaPP	Robust PCA by Projection Pursuit
PCDSpline	Semiparametric regression analysis of panel count data using monotone splines
pcev	Principal Component of Explained Variance
pcg	Preconditioned Conjugate Gradient Algorithm for solving Ax=b
PCGSE	Principal Component Gene Set Enrichment
pch	Piecewise Constant Hazards Models for Censored and Truncated Data
PCICt	Implementation of POSIXct work-alike for 365 and 360 day calendars
pcIRT	IRT Models for Polytomous and Continuous Item Responses
PCIT	Partial Correlation Coefficient with Information Theory
pcnetmeta	Patient-Centered Network Meta-Analysis
pco	Panel Cointegration Tests
PCovR	Principal Covariates Regression
PCPS	Principal Coordinates of Phylogenetic Structure
PCS	Calculate the probability of correct selection (PCS)
pcse	Panel-Corrected Standard Error Estimation in R
pdc	Permutation Distribution Clustering
pdfCluster	Cluster analysis via nonparametric density estimation
pdfetch	Fetch Economic and Financial Time Series Data from Public Sources
pdftables	Programmatic Conversion of PDF Tables
pdftools	Extract Text and Data from PDF Documents
pdist	Partitioned Distance Function
pdmod	Proximal/distal modeling framework for Pavlovian conditioning phenomena
pdolsms	Panel Dynamic OLS Estimation of Cointegrating Vectors
PDQutils	PDQ Functions via Gram Charlier, Edgeworth, and Cornish Fisher Approximations
pdR	Panel Data Regression: Threshold Model and Unit Root Tests
PDSCE	Positive definite sparse covariance estimators
peacots	Periodogram Peaks in Correlated Time Series
peakPick	Peak Picking Methods Inspired by Biological Data
Peaks	Peaks
pear	Package for Periodic Autoregression Analysis
pearson7	Maximum Likelihood Inference for the Pearson VII Distribution with Shape Parameter 3/2
PearsonDS	Pearson Distribution System
PearsonICA	Independent component analysis using score functions from the Pearson system
pec	Prediction Error Curves for Risk Prediction Models in Survival Analysis
pedantics	Functions to facilitate power and sensitivity analyses for genetic studies of natural populations
PedCNV	An implementation for association analysis with CNV data
pedgene	Gene-Level Statistics for Pedigree Data
pedigree	Pedigree functions
pedigreemm	Pedigree-based mixed-effects models
pedometrics	Pedometric Tools and Techniques
pegas	Population and Evolutionary Genetics Analysis System
PEIP	Geophysical Inverse Theory and Optimization
PEMM	A Penalized EM algorithm incorporating missing-data mechanism
pems.utils	Portable Emissions (and Other Mobile) Measurement System Utilities
penalized	L1 (lasso and fused lasso) and L2 (ridge) penalized estimation in GLMs and in the Cox model
penalizedLDA	Penalized Classification using Fisher's Linear Discriminant
penalizedSVM	Feature Selection SVM using penalty functions
pencopula	Flexible Copula Density Estimation with Penalized Hierarchical B-Splines
PenCoxFrail	Regularization in Cox Frailty Models
pendensity	Density Estimation with a Penalized Mixture Approach
penDvine	Flexible Pair-Copula Estimation in D-Vines using Bivariate Penalized Splines
penMSM	Estimating Regularized Multi-state Models Using L1 Penalties
pensim	Simulation of high-dimensional data and parallelized repeated penalized regression
peperr	Parallelised Estimation of Prediction Error
peplib	Peptide Library Analysis Methods
PepPrep	Insilico peptide mutation, digestion and homologous comparison
peptider	Evaluation of Diversity in Nucleotide Libraries
Peptides	Calculate Indices and Theoretical Properties of Protein Sequences
pequod	Moderated Regression Package
perARMA	Periodic Time Series Analysis
Perc	Using Percolation and Conductance to Find Information Flow Certainty in a Direct Network
PerFit	Person Fit
PerfMeas	PerfMeas: Performance Measures for ranking and classification tasks
PerformanceAnalytics	Econometric tools for performance and risk analysis
performanceEstimation	An Infra-Structure for Performance Estimation of Predictive Models
perm	Exact or Asymptotic permutation tests
PermAlgo	Permutational Algorithm to Simulate Survival Data
PerMallows	Permutations and Mallows Distributions
permGPU	Using GPUs in Statistical Genomics
permPATH	Permutation-Based Gene Expression Pathway Analyses
permubiome	A Permutation Based Test for Biomarker Discovery in Microbiome Data
permutations	Permutations of a Finite Set
permute	Functions for Generating Restricted Permutations of Data
perry	Resampling-based prediction error estimation for regression models
persiandictionary	English to Persian dictionary
personograph	Pictographic Representation of Treatment Effects
perspectev	Permutation of Species During Turnover Events
perturb	Tools for evaluating collinearity
pesticides	Analysis of single serving and composite pesticide residue measurements
PET	Simulation and Reconstruction of PET Images
pez	Phylogenetics for the Environmental Sciences
pgam	Poisson-Gamma Additive Models
PGEE	Penalized Generalized Estimating Equations in High-Dimension
PGICA	Parallel Group ICA Algorithm
pgirmess	Data Analysis in Ecology
pglm	panel generalized linear model
pGLS	Generalized Least Square in comparative Phylogenetics
PGM2	Recursive method for construction of nested resolvable designs and uniform designs associated
pgmm	Parsimonious Gaussian Mixture Models
pgnorm	The p-Generalized Normal Distribution
PGRdup	Discover Probable Duplicates in Plant Genetic Resources Collections
pgs	Precision of Geometric Sampling
ph2bayes	Bayesian Single-Arm Phase II Designs
phalen	Phalen Algorithms and Functions
phangorn	Phylogenetic Analysis in R
PharmacoGx	Analysis of Large-Scale Pharmacogenomic Data
PharmPow	Pharmacometric Power calculations for mixed study designs
phaseR	Phase Plane Analysis of One and Two Dimensional Autonomous ODE Systems
PhaseType	Inference for Phase-type Distributions
phcfM	Modelling anthropogenic deforestation
pheatmap	Pretty Heatmaps
phenability	Nonparametric Stability Analysis
phenex	Auxiliary Functions for Phenological Data Analysis
PHENIX	Phenotypic Integration Index
phenmod	Auxiliary functions for phenological data processing, modelling and result handling
pheno	Auxiliary functions for phenological data analysis
pheno2geno	High-Throughput Generation of Genetic Markers and Maps from Molecular Phenotypes for Crosses Between Inbred Strains
phenology	Tools to Manage a Parametric Function that Describes Phenology
PHeval	Evaluation of the Proportional Hazards Assumption with a Standardized Score Process
phia	Post-Hoc Interaction Analysis
phmm	Proportional Hazards Mixed-effects Model (PHMM)
phonenumber	Convert Letters to Numbers and Back as on a Telephone Keypad
phonics	Phonetic Spelling Algorithms
phonR	Tools for Phoneticians and Phonologists
phonTools	Tools for Phonetic and Acoustic Analyses
photobiology	Photobiological Calculations
photobiologyWavebands	Waveband Definitions for UV, VIS, and IR Radiation
phreeqc	R Interface to Geochemical Modeling Software
phtt	Panel Data Analysis with Heterogeneous Time Trends
PhViD	PhViD: an R package for PharmacoVigilance signal Detection
Phxnlme	Run Phoenix NLME and Perform Post-Processing
PhyActBedRest	Marks periods of sleep in Actigraph accelerometer data
phyclust	Phylogenetic Clustering (Phyloclustering)
phyext2	An Extension (for Package 'SigTree') of Some of the Classes in Package 'phylobase'
phylin	Spatial Interpolation of Genetic Data
phylobase	Base Package for Phylogenetic Structures and Comparative Data
phyloclim	Integrating Phylogenetics and Climatic Niche Modeling
phylocurve	Phylogenetic Comparative Methods for High-Dimensional Traits
PHYLOGR	Functions for Phylogenetically Based Statistical Analyses
phyloland	Modelling Competitive Exclusion and Limited Dispersal in a Statistical Phylogeographic Framework
phylolm	Phylogenetic Linear Regression
PhyloMeasures	Fast and Exact Algorithms for Computing Phylogenetic Biodiversity Measures
phylometrics	Estimating Statistical Errors of Phylogenetic Metrics
phylosignal	Exploring the Phylogenetic Signal in Continuous Traits
phylotools	Phylogenetic tools for Eco-phylogenetics
phyndr	Matches Tip and Trait Data
phyreg	Implements the Phylogenetic Regression of Grafen (1989)
PhysicalActivity	Process Physical Activity Accelerometer Data
physiology	Calculate Physiological Characteristics of Adults and Children
PhySortR	A Fast, Flexible Tool for Sorting Phylogenetic Trees
phytools	Phylogenetic Tools for Comparative Biology (and Other Things)
phytotools	Phytoplankton Production Tools
pi0	Estimating the Proportion of True Null Hypotheses for FDR
picante	R tools for integrating phylogenies and ecology
picasso	Pathwise Calibrated Sparse Shooting Algorithm
pid	Process Improvement using Data
piecewiseSEM	Piecewise Structural Equation Modeling
PIGE	Self contained gene set analysis for gene- and pathway-environment interaction analysis
PIGShift	Polygenic Inverse Gamma Shifts
Pijavski	Global Univariate Minimization
pinfsc50	Sequence ("Fasta"), Annotation ("Gff") and Variants ("Vcf") for 17 Samples of "P. Infestans" and 1 "P. Mirabilis"
pingr	Check If a Remote Computer is Up
pinnacle.API	A Wrapper for the Pinnacle Sports API
pipe.design	Dual-Agent Dose Escalation for Phase I Trials using the PIPE Design
pipeR	Multi-Paradigm Pipeline Implementation
PIPS	Predicted Interval Plots
pitchRx	Tools for Harnessing 'MLBAM' 'Gameday' Data and Visualizing 'pitchfx'
PivotalR	A Fast, Easy-to-use Tool for Manipulating Tables in Databases and A Wrapper of MADlib
pixiedust	Tables so Beautifully Fine-Tuned You Will Believe It's Magic
pixmap	Bitmap Images (“Pixel Maps”)
PK	Basic Non-Compartmental Pharmacokinetics
pkgconfig	Private Configuration for 'R' Packages
pkgKitten	Create Simple Packages Which Do not Upset R Package Checks
pkgmaker	Package development utilities
PKgraph	Model diagnostics for population pharmacokinetic models
PKI	Public Key Infrastucture for R Based on the X.509 Standard
PKNCA	Perform Pharmacokinetic Non-Compartmental Analysis
PKPDmodels	Pharmacokinetic/pharmacodynamic models
PKreport	A reporting pipeline for checking population pharmacokinetic model assumption
pks	Probabilistic Knowledge Structures
pla	Parallel Line Assays
plan	Tools for project planning
planar	Multilayer Optics
planor	Generation of Regular Factorial Designs
plantecophys	Modelling and Analysis of Leaf Gas Exchange Data
plaqr	Partially Linear Additive Quantile Regression
PlayerRatings	Dynamic Updating Methods for Player Ratings Estimation
playwith	A GUI for interactive plots using GTK+
plfm	Probabilistic Latent Feature Analysis
plfMA	A GUI to View, Design and Export Various Graphs of Data
plgp	Particle Learning of Gaussian Processes
PLIS	Multiplicity control using Pooled LIS statistic
plm	Linear Models for Panel Data
plmDE	Additive partially linear models for differential gene expression analysis
plmm	Partially Linear Mixed Effects Model
pln	Polytomous logit-normit (graded logistic) model estimation
PLordprob	Multivariate Ordered Probit Model via Pairwise Likelihood
plot3D	Plotting Multi-Dimensional Data
plot3Drgl	Plotting Multi-Dimensional Data - Using 'rgl'
plotGoogleMaps	Plot Spatial or Spatio-Temporal Data Over Google Maps
plotKML	Visualization of Spatial and Spatio-Temporal Objects in Google Earth
plotly	Create Interactive Web Graphics via 'plotly.js'
plotMCMC	MCMC Diagnostic Plots
plotmo	Plot a Model's Response and Residuals
plotpc	Plot Principal Component Histograms Around a Scatter Plot
PlotPrjNetworks	Useful Networking Tools for Project Management
PlotRegionHighlighter	Creates an envelope that surrounds a set of points plotted in a two dimensional space
plotrix	Various Plotting Functions
plotROC	Generate Useful ROC Curve Charts for Print and Interactive Use
plotSEMM	Graphing Nonlinear Relations Among Latent Variables from Structural Equation Mixture Models
plRasch	Log Linear by Linear Association models and Rasch family models by pseudolikelihood estimation
PLRModels	Statistical inference in partial linear regression models
pls	Partial Least Squares and Principal Component Regression
PLSbiplot1	The Partial Least Squares (PLS) Biplot
plsdepot	Partial Least Squares (PLS) Data Analysis Methods
plsdof	Degrees of Freedom and Statistical Inference for Partial Least Squares Regression
plsgenomics	PLS Analyses for Genomics
plspm	Tools for Partial Least Squares Path Modeling (PLS-PM)
plspm.formula	Formula Based PLS Path Modeling
plsRbeta	Partial Least Squares Regression for Beta Regression Models
plsRcox	Partial Least Squares Regression for Cox Models and Related Techniques
plsRglm	Partial Least Squares Regression for Generalized Linear Models
plugdensity	Plug-in Kernel Density Estimation
plumbr	Mutable and dynamic data models
plus	Penalized Linear Unbiased Selection
plusser	A Google+ Interface for R
plyr	Tools for Splitting, Applying and Combining Data
PMA	Penalized Multivariate Analysis
pmc	Phylogenetic Monte Carlo
pmcgd	pmcgd
pmclust	Parallel Model-Based Clustering using Expectation-Gathering-Maximization Algorithm for Finite Mixture Gaussian Model
PMCMR	Calculate Pairwise Multiple Comparisons of Mean Rank Sums
pmg	Poor Man's GUI
pmhtutorial	Minimal Working Examples for Particle Metropolis-Hastings
pmlr	Penalized Multinomial Logistic Regression
pmml	Generate PMML for Various Models
pmmlTransformations	Transforms Input Data from a PMML Perspective
pmr	Probability Models for Ranking Data
pnea	Parametric Network Enrichment Analysis
pnf	Prime Numbers and Integer Factorization
png	Read and write PNG images
pnmtrem	Probit-Normal Marginalized Transition Random Effects Models
pnn	Probabilistic neural networks
pocrm	Dose Finding in Drug Combination Phase I Trials Using PO-CRM
pogit	Bayesian Variable Selection for a Poisson-Logistic Model
PogromcyDanych	PogromcyDanych / DataCrunchers is the Masive Online Open Course that Brings R and Statistics to the People
poibin	The Poisson Binomial Distribution
PoiClaClu	Classification and clustering of sequencing data based on a Poisson model
poilog	Poisson lognormal and bivariate Poisson lognormal distribution
pointdensityP	Point Density for Geospatial Data
pointRes	Analyzing Pointer Years and Components of Resilience
PoisBinNonNor	Data Generation with Poisson, Binary and Continuous Components
PoisBinOrd	Data Generation with Poisson, Binary and Ordinal Components
PoisBinOrdNonNor	Generation of up to Four Different Types of Variables
PoisBinOrdNor	Data Generation with Poisson, Binary, Ordinal and Normal Components
poisDoubleSamp	Confidence Intervals with Poisson Double Sampling
PoisNonNor	Simultaneous Generation of Count and Continuous Data
PoisNor	Simultaneous generation of multivariate data with Poisson and normal marginals
poisson	Simulating Homogenous & Non-Homogenous Poisson Processes
poisson.glm.mix	Fit high dimensional mixtures of Poisson GLMs
PoissonSeq	Significance analysis of sequencing data based on a Poisson log linear model
poistweedie	Poisson-Tweedie exponential family models
poLCA	Polytomous variable Latent Class Analysis
polidata	Political Data Interface in R
pollstR	Client for the HuffPost Pollster API
polspline	Polynomial Spline Routines
polyaAeppli	Implementation of the Polya-Aeppli distribution
polyapost	Simulating from the Polya Posterior
polychaosbasics	Sensitivity Indexes Calculated from Polynomial Chaos Expansions
polyclip	Polygon Clipping
polycor	Polychoric and Polyserial Correlations
polyCub	Cubature over Polygonal Domains
polyfreqs	Bayesian Population Genomics in Autopolyploids
polynom	A Collection of Functions to Implement a Class for Univariate Polynomial Manipulations
PolynomF	Polynomials in R
PolyPatEx	Paternity Exclusion in Autopolyploid Species
polysat	Tools for Polyploid Microsatellite Analysis
polySegratio	Simulate and test marker dosage for dominant markers in autopolyploids
polySegratioMM	Bayesian mixture models for marker dosage in autopolyploids
polywog	Bootstrapped Basis Regression with Oracle Model Selection
pom	POM - Patch Occupancy Models
Pomic	Pomic
pomp	Statistical Inference for Partially Observed Markov Processes
pooh	Partial Orders and Relations
popbio	Construction and Analysis of Matrix Population Models
popdemo	Provides Tools For Demographic Modelling Using Projection Matrices
PopED	Population (and Individual) Optimal Experimental Design
popEpi	Functions for Epidemiological Analysis using Population Data
PopGenKit	Useful functions for (batch) file conversion and data resampling in microsatellite datasets
PopGenome	An Efficient Swiss Army Knife for Population Genomic Analyses
PopGenReport	A Simple Framework to Analyse Population Genetic Data
popgraph	This is an R package that constructs and manipulates population graphs
popKorn	For interval estimation of mean of selected populations
poplite	Tools for Simplifying the Population and Querying of SQLite Databases
poppr	Genetic Analysis of Populations with Mixed Reproduction
popprxl	Read GenAlEx Files Directly from Excel
popRange	popRange: A spatially and temporally explicit forward genetic simulator
popReconstruct	Reconstruct Human Populations of the Recent Past
popsom	Routines for Constructing and Evaluating Self-Organizing Maps
population	Models for Simulating Populations
PopVar	Genomic Breeding Tools: Genetic Variance Prediction and Cross-Validation
portes	Portmanteau Tests for Univariate and Multivariate Time Series Models
portfolio	Analysing equity portfolios
PortfolioAnalytics	Portfolio Analysis, Including Numerical Methods for Optimization of Portfolios
PortfolioEffectEstim	High Frequency Price Estimators by PortfolioEffect
PortfolioEffectHFT	High Frequency Portfolio Analytics by PortfolioEffect
portfolioSim	Framework for simulating equity portfolio strategies
PortRisk	Portfolio Risk Analysis
potts	Markov Chain Monte Carlo for Potts Models
PottsUtils	Utility Functions of the Potts Models
powell	Powell's UObyQA algorithm
PoweR	Computation of Power and Level Tables for Hypothesis Tests
Power2Stage	Power and Sample-Size Distribution of 2-Stage Bioequivalence Studies
powerAnalysis	Power analysis in experimental design
powerGWASinteraction	Power Calculations for GxE and GxG Interactions for GWAS
poweRlaw	Analysis of Heavy Tailed Distributions
powerMediation	Power/Sample Size Calculation for Mediation Analysis
powerpkg	Power analyses for the affected sib pair and the TDT design
powerplus	Exponentiation Operations
powerSurvEpi	Power and Sample Size Calculation for Survival Analysis of Epidemiological Studies
PowerTOST	Power and Sample Size Based on Two One-Sided t-Tests (TOST) for (Bio)Equivalence Studies
PP	Person Parameter estimation
ppcor	Partial and Semi-Partial (Part) Correlation
ppiPre	Predict Protein-Protein Interactions Based on Functional and Topological Similarities
ppls	Penalized Partial Least Squares
ppmlasso	Point Process Models with LASSO Penalties
pps	Functions for PPS sampling
PPtree	Projection pursuit classification tree
PPtreeViz	Projection Pursuit Classification Tree Visualization
pqantimalarials	web tool for estimating under-five deaths caused by poor-quality antimalarials in sub-Saharan Africa
prabclus	Functions for Clustering of Presence-Absence, Abundance and Multilocus Genetic Data
pracma	Practical Numerical Math Functions
PracTools	Tools for Designing and Weighting Survey Samples
pragma	Provides a pragma / directive / keyword syntax for R
prais	Prais-Winsten Estimation Procedure for AR(1) Serial Correlation
praise	Praise Users
praktikum	Kvantitatiivsete meetodite praktikumi asjad / Functions used in the course "Quantitative methods in behavioural sciences" (SHPH.00.004), University of Tartu
prc	Paired Response Curve
prcbench	Testing Workbench for Precision-Recall Curves
prclust	Penalized Regression-Based Clustering Method
precintcon	Precipitation Intensity, Concentration and Anomaly Analysis
precrec	Calculate Accurate Precision-Recall and Receiver Operator Characteristics Curves
PredictABEL	Assessment of Risk Prediction Models
PredictiveRegression	Prediction Intervals for Three Basic Statistical Models
predictmeans	Calculate Predicted Means for Linear Models
predmixcor	Classification rule based on Bayesian mixture models with feature selection bias corrected
prefmod	Utilities to Fit Paired Comparison Models for Preferences
PReMiuM	Dirichlet Process Bayesian Clustering, Profile Regression
prepdat	Preparing Experimental Data for Statistical Analysis
preprocomb	Tools for Preprocessing Combinations
preproviz	Tools for Visualization of Interdependent Data Quality Issues
prereg	R Markdown Template to Preregister Scientific Studies
PresenceAbsence	Presence-Absence Model Evaluation
presens	R Interface for PreSens Fiber Optic Data
preseqR	Predicting Species Accumulation Curves
prettyGraphs	publication-quality graphics
prettymapr	Scale Bar, North Arrow, and Pretty Margins in R
prettyR	Pretty Descriptive Stats
prettyunits	Pretty, Human Readable Formatting of Quantities
prevalence	Tools for Prevalence Assessment Studies
PrevMap	Geostatistical Modelling of Spatially Referenced Prevalence Data
prevR	Estimating Regional Trends of a Prevalence from a DHS
pRF	Permutation Significance for Random Forests
prim	Patient Rule Induction Method (PRIM)
primer	Functions and data for A Primer of Ecology with R
primerTree	Visually Assessing the Specificity and Informativeness of Primer Pairs
primes	Generate and Test for Prime Numbers
PRIMsrc	PRIM Survival Regression Classification
princurve	Fits a Principal Curve in Arbitrary Dimension
prinsimp	Finding and plotting simple basis vectors for multivariate data
prism	Access Data from the Oregon State Prism Climate Project
PRISMA	Protocol Inspection and State Machine Analysis
PrivateLR	Differentially Private Regularized Logistic Regression
prLogistic	Estimation of Prevalence Ratios using Logistic Models
pro	Point-Process Response Model for Optogenetics
prob	Elementary Probability on Finite Sample Spaces
probemod	Statistical Tools for Probing Moderation Effects
probFDA	Probabilistic Fisher Discriminant Analysis
ProbForecastGOP	Probabilistic weather forecast using the GOP method
ProbitSpatial	Probit with Spatial Dependence, SAR and SEM Models
probsvm	probsvm: Class probability estimation for Support Vector Machines
ProbYX	Inference for the Stress-Strength Model R = P(Y<X)
pROC	Display and Analyze ROC Curves
ProDenICA	Product Density Estimation for ICA using tilted Gaussian density estimates
prodlim	Product-Limit Estimation for Censored Event History Analysis
PROFANCY	The package can prioritize candidate disease metabolites based on global functional relationships between metabolites in the context of metabolic pathways
profdpm	Profile Dirichlet Process Mixtures
ProfessR	Grades Setting and Exam Maker
ProfileLikelihood	Profile Likelihood for a Parameter in Commonly Used Statistical Models
profileModel	Tools for profiling inference functions for various model classes
profileR	Profile Analysis of Multivariate Data in R
profilr	Quickly Profile Data in R
profr	An alternative display for profiling information
proftools	Profile Output Processing Tools for R
progenyClust	Finding the Optimal Cluster Number Using Progeny Clustering
ProgGUIinR	support package for "Programming Graphical User Interfaces in R"
prognosticROC	Prognostic ROC curves for evaluating the predictive capacity of a binary test
progress	Terminal Progress Bars
proj4	A simple interface to the PROJ.4 cartographic projections library
ProjectTemplate	Automates the creation of new statistical analysis projects
ProNet	Biological Network Construction, Visualization and Analyses
propagate	Propagation of Uncertainty
PropCIs	Various confidence interval methods for proportions
PropClust	Propensity Clustering and Decomposition
prop.comb.RR	Analyzing Combination of Proportions and Relative Risk
properties	Parse Java Properties Files for R Service Bus Applications
proportion	Inference on Single Binomial Proportion and Bayesian Computations
propOverlap	Feature (gene) selection based on the Proportional Overlapping Scores
PropScrRand	Propensity score methods for assigning treatment in randomized trials
prospectr	Miscellaneous functions for processing and sample selection of vis-NIR diffuse reflectance data
ProteinDescriptors	Generates Various Protein Descriptors for Machine Learning Algorithms
proteomicdesign	Optimization of a multi-stage proteomic study
proteomics	Statistical Analysis of High Throughput Proteomics Data
protiq	Protein (identification and) quantification based on peptide evidence
proto	Prototype object-based programming
protoclass	Interpretable classification with prototypes
protoclust	Hierarchical Clustering with Prototypes
PROTOLIDAR	PRocess TOol LIdar DAta in R
proton	The Proton Game
prototest	Inference on Prototypes from Clusters of Features
protr	Generating Various Numerical Representation Schemes of Protein Sequence
ProTrackR	Manipulate and Play 'ProTracker' Modules
protViz	Visualizing and Analyzing Mass Spectrometry Related Data in Proteomics
provenance	Statistical Toolbox for Sedimentary Provenance Analysis
proxy	Distance and Similarity Measures
prozor	Minimal Protein Set Explaining Peptide Spectrum Matches
PRROC	Precision-Recall and ROC Curves for Weighted and Unweighted Data
pRSR	Test of periodicity using response surface regression
pryr	Tools for Computing on the Language
PSAboot	Bootstrapping for Propensity Score Analysis
PSAgraphics	Propensity Score Analysis Graphics
psbcGroup	Penalized Parametric and Semiparametric Bayesian Survival Models with Shrinkage and Grouping Priors
PSCBS	Analysis of Parent-Specific DNA Copy Numbers
pscl	Political Science Computational Laboratory, Stanford University
pscore	Standardizing Physiological Composite Risk Endpoints
psd	Adaptive, Sine-Multitaper Power Spectral Density Estimation
psData	Download Regularly Maintained Political Science Data Sets
pse	Parameter Space Exploration with Latin Hypercubes
pseudo	Pseudo - observations
pseval	Methods for Evaluating Principal Surrogates of Treatment Response
PSF	Algorithm for Pattern Sequence Based Forecasting
psgp	Projected Spatial Gaussian Process (psgp) methods
pSI	Specificity Index Statistic
psidR	Build Panel Data Sets from PSID Raw Data
PsiHat	Several Local False Discovery Rate Estimators
PSM	Non-Linear Mixed-Effects modelling using Stochastic Differential Equations
pso	Particle Swarm Optimization
psoptim	Particle Swarm Optimization
pspearman	Spearman's rank correlation test
pspline	Penalized Smoothing Splines
pssm	Piecewise Exponential Model for Time to Progression and Time from Progression to Death
PST	Probabilistic Suffix Trees and Variable Length Markov Chains
PsumtSim	Simulations of grouped responses relative to baseline
psy	Various procedures used in psychometry
psych	Procedures for Psychological, Psychometric, and Personality Research
psychometric	Applied Psychometric Theory
psychomix	Psychometric Mixture Models
psychotools	Infrastructure for Psychometric Modeling
psychotree	Recursive Partitioning Based on Psychometric Models
psyphy	Functions for analyzing psychophysical data in R
psytabs	Produce well-formatted tables for psychological research
PTAk	Principal Tensor Analysis on k Modes
PTE	Personalized Treatment Evaluator
ptinpoly	Point-In-Polyhedron Test (2D and 3D)
PtProcess	Time Dependent Point Process Modelling
ptw	Parametric Time Warping
ptycho	Bayesian Variable Selection with Hierarchical Priors
PubBias	Performs simulation study to look for publication bias, using a technique described by Ioannidis and Trikalinos; Clin Trials. 2007;4(3):245-53
pubmed.mineR	Text Mining of PubMed Abstracts
PubMedWordcloud	PubMed Word Clouds
pubprint	Printing Results of Statistical Computing in a Publishable Way
pullword	R Interface to Pullword Service
pumilioR	Pumilio in R
PurBayes	Bayesian Estimation of Tumor Purity and Clonality
purge	Purge Training Data from Models
purrr	Functional Programming Tools
pushoverr	Send push notifications using Pushover
PVAClone	Population Viability Analysis with Data Cloning
pvar	Calculation and Application of p-variation
pvclass	P-Values for Classification
pvclust	Hierarchical Clustering with P-Values via Multiscale Bootstrap Resampling
PVR	Computes phylogenetic eigenvectors regression (PVR) and phylogenetic signal-representation curve (PSR) (with null and Brownian expectations)
pvrank	Rank Correlations
pvsR	An R package to interact with the Project Vote Smart API for scientific research
PWD	Time Series Regression Using the Power Weighted Densities (PWD) Approach
pweight	P-Value Weighting
pwr	Basic Functions for Power Analysis
PwrGSD	Power in a Group Sequential Design
pwrRasch	Statistical Power Simulation for Testing the Rasch Model
pwt	Penn World Table (Versions 5.6, 6.x, 7.x)
pwt8	Penn World Table (Version 8.x)
pxR	PC-Axis with R
pxweb	R Interface to the PX-Web/PC-Axis API
pycno	Pycnophylactic Interpolation
pyramid	Functions to draw population pyramid
pystr	Python String Methods in R
PythonInR	Use Python from Within R
0	
qap	Heuristics for the Quadratic Assignment Problem (QAP)
QCA	Qualitative Comparative Analysis
QCA3	Yet another package for Qualitative Comparative Analysis
QCAfalsePositive	Tests for Type I Error in Qualitative Comparative Analysis (QCA)
QCAGUI	Modern Functions for Qualitative Comparative Analysis
QCApro	Professional Functionality for Performing and Evaluating Qualitative Comparative Analysis
QCAtools	Helper Functions for QCA in R
qcc	Quality Control Charts
QCGWAS	Quality Control of Genome Wide Association Study results
qclust	Robust Estimation of Gaussian Mixture Models
qcr	Quality control and reliability
QCSIS	Sure Independence Screening via Quantile Correlation and Composite Quantile Correlation
qdap	Bridging the Gap Between Qualitative Data and Quantitative Analysis
qdapDictionaries	Dictionaries and Word Lists for the 'qdap' Package
qdapRegex	Regular Expression Removal, Extraction, and Replacement Tools
qdapTools	Tools for the 'qdap' Package
qdm	Fitting a Quadrilateral Dissimilarity Model to Same-Different Judgments
QFRM	Pricing of Vanilla and Exotic Option Contracts
qgraph	Graph Plotting Methods, Psychometric Data Visualization and Graphical Model Estimation
qgtools	Tools for Quantitative Genetics Data Analyses
QICD	Estimate the Coefficients for Non-Convex Penalized Quantile Regression Model by using QICD Algorithm
qicharts	Quality Improvement Charts
qiimer	Work with QIIME Output Files in R
qlcData	Processing Data for Quantitative Language Comparison (QLC)
qlcMatrix	Utility Sparse Matrix Functions for Quantitative Language Comparison
qlcVisualize	Visualization for Quantitative Language Comparison (QLC)
qLearn	Estimation and inference for Q-learning
qmap	Statistical transformations for post-processing climate model output
qmethod	Analysis of Subjective Perspectives Using Q Methodology
qmrparser	Parser combinator in R
QoLR	Analysis of Health-Related Quality of Life in Oncology
qpcR	Modelling and analysis of real-time PCR data
qPCR.CT	qPCR data analysis and plot package
QPot	Quasi-Potential Analysis for Stochastic Differential Equations
qqman	Q-Q and manhattan plots for GWAS data
qqtest	Self Calibrating Quantile-Quantile Plots for Visual Testing
qrage	Tools that Create D3 JavaScript Force Directed Graph from R
qrcm	Quantile Regression Coefficients Modeling
qrcode	QRcode Generator for R
qrfactor	Simultaneous simulation of Q and R mode factor analyses with Spatial data
qrjoint	Joint Estimation in Linear Quantile Regression
qrLMM	Quantile Regression for Linear Mixed-Effects Models
QRM	Provides R-Language Code to Examine Quantitative Risk Management Concepts
qrmtools	Tools for Quantitative Risk Management
qrng	(Randomized) Quasi-Random Number Generators
qrNLMM	Quantile Regression for Nonlinear Mixed-Effects Models
qrnn	Quantile Regression Neural Network
QSARdata	Quantitative Structure Activity Relationship (QSAR) Data Sets
qtbase	Interface Between R and Qt
qte	Quantile Treatment Effects
qtl	Tools for Analyzing QTL Experiments
qtlbook	Datasets for the R/qtl Book
qtlc	Densitometric Analysis of Thin-Layer Chromatography Plates
qtlcharts	Interactive Graphics for QTL Experiments
qtlDesign	Design of QTL experiments
qtlhot	Inference for QTL Hotspots
qtlmt	Tools for Mapping Multiple Complex Traits
qtlnet	Causal Inference of QTL Networks
QTLRel	Tools for Mapping of Quantitative Traits of Genetically Related Individuals and Calculating Identity Coefficients from Pedigrees
Qtools	Utilities for Quantiles
qtpaint	Qt-Based Painting Infrastructure
qtutils	Miscellaneous Qt-based utilities
QuACN	QuACN: Quantitative Analysis of Complex Networks
quad	Exact permutation moments of quadratic form statistics
quadprog	Functions to solve Quadratic Programming Problems
quadrupen	Sparsity by Worst-Case Quadratic Penalties
qualCI	Causal Inference with Qualitative and Ordinal Information on Outcomes
QualInt	Test for Qualitative Interactions
qualityTools	Statistical Methods for Quality Science
qualV	Qualitative Validation Methods
qualvar	Implements Indices of Qualitative Variation Proposed by Wilcox (1973)
Quandl	API Wrapper for Quandl.com
quantable	Streamline Descriptive Analysis of Quantitative Data Matrices
quantchem	Quantitative chemical analysis: calibration and evaluation of results
quanteda	Quantitative Analysis of Textual Data
quantification	Quantification of Qualitative Survey Data
QuantifQuantile	Estimation of Conditional Quantiles using Optimal Quantization
quantileDA	Quantile Classifier
quantmod	Quantitative Financial Modelling Framework
QuantPsyc	Quantitative Psychology Tools
quantreg	Quantile Regression
quantregForest	Quantile Regression Forests
quantregGrowth	Growth Charts via Regression Quantiles
quantspec	Quantile-Based Spectral Analysis of Time Series
QuantumClone	Clustering Mutations using High Throughput Sequencing (HTS) Data
QuasiSeq	Analyzing RNA Sequencing Count Tables Using Quasi-Likelihood
questionr	Functions to Make Surveys Processing Easier
queueing	Analysis of Queueing Networks and Models
QUIC	Regularized sparse inverse covariance matrix estimation
quickmapr	Quickly Map and Explore Spatial Data
quickpsy	Fits Psychometric Functions for Multiple Groups
quickReg	Build Regression Models Quickly and Display the Results Using 'ggplot2'
quint	Qualitative Interaction Trees
quipu	Summary charts of micro satellite profiles for a set of biological samples
Quor	Quantile Ordering
qut	Quantile Universal Threshold
qVarSel	Variables Selection for Clustering and Classification
qvcalc	Quasi Variances for Factor Effects in Statistical Models
qwraps2	Quick Wraps 2
QZ	Generalized Eigenvalues and QZ Decomposition
0	
R0	Estimation of R0 and Real-Time Reproduction Number from Epidemics
R1magic	Compressive Sampling: Sparse Signal Recovery Utilities
R2admb	'ADMB' to R Interface Functions
R2BayesX	Estimate Structured Additive Regression Models with BayesX
R2Cuba	Multidimensional Numerical Integration
r2d2	Bivariate (Two-Dimensional) Confidence Region and Frequency Distribution
r2dRue	2d Rain Use Efficience model
R2G2	Converting R CRAN outputs into Google Earth
R2GUESS	Wrapper Functions for GUESS
R2HTML	HTML exportation for R objects
R2jags	Using R to Run 'JAGS'
r2lh	R to LaTeX and HTML
R2MLwiN	Running 'MLwiN' from Within R
R2OpenBUGS	Running OpenBUGS from R
R2PPT	Simple R Interface to Microsoft PowerPoint using rcom or RDCOMClient
R2STATS	A GTK GUI for fitting and comparing GLM and GLMM in R
r2stl	r2stl, R package for visualizing data using a 3D printer
R2SWF	Convert R Graphics to Flash Animations
R2wd	Write MS-Word documents from R
R2WinBUGS	Running 'WinBUGS' and 'OpenBUGS' from 'R' / 'S-PLUS'
R330	An R package for Stats 330
R4CouchDB	An R Convenience Layer for CouchDB
R4dfp	4dfp MRI Image Read and Write Routines
r4ss	R Code for Stock Synthesis
R6	Classes with Reference Semantics
race	Racing methods for the selection of the best
RAD	Fit RAD models to biological data
RADami	R Package for Phylogenetic Analysis of RADseq Data
RADanalysis	Normalization and Study of Rank Abundance Distributions
radar	Fundamental Formulas for Radar
radarchart	Radar Chart from Chart.js
radiomics	Radiomic Image Processing Toolbox
RadioSonde	Tools for plotting skew-T diagrams and wind profiles
radir	Inverse-Regression Estimation of Radioactive Doses
RadOnc	Analytical Tools for Radiation Oncology
RadTran	Radon and Soil Gas Transport in 2D Porous Medium
RAdwords	Loading Google Adwords Data into R
rafalib	Convenience Functions for Routine Data Exploration
rags2ridges	Ridge Estimation of Precision Matrices from High-Dimensional Data
RAHRS	Data Fusion Filters for Attitude Heading Reference System (AHRS) with Several Variants of the Kalman Filter and the Mahoney and Madgwick Filters
rainbow	Rainbow Plots, Bagplots and Boxplots for Functional Data
raincpc	Obtain and Analyze Rainfall Data from the Climate Prediction Center
rainfreq	Rainfall Frequency (Design Storm) Estimates from the US National Weather Service
rAltmetric	Retrieves Altmerics Data For Any Published Paper From Altmetric.com
RAM	R for Amplicon-Sequencing-Based Microbial-Ecology
Rambo	The Random Subgraph Model
rAmCharts	JavaScript Charts API Tool
Ramd	Tools For Managing File/function Dependencies In R
ramify	Additional Matrix Functionality
RAMP	Regularized Generalized Linear Models with Interaction Effects
RAMpath	Structural Equation Modeling using the reticular action model (RAM) Notation
ramps	Bayesian Geostatistical Modeling with RAMPS
ramsvm	Reinforced Angle-Based Multicategory Support Vector Machines
randaes	Random number generator based on AES cipher
randNames	Package Provides Access to Fake User Data
random	True Random Numbers using RANDOM.ORG
randomcoloR	Generate Attractive Random Colors
RandomFields	Simulation and Analysis of Random Fields
RandomFieldsUtils	Utilities for the Simulation and Analysis of Random Fields
randomForest	Breiman and Cutler's Random Forests for Classification and Regression
randomForest.ddR	Distributed 'randomForest' for Big Data using 'ddR' API
randomForestSRC	Random Forests for Survival, Regression and Classification (RF-SRC)
randomGLM	Random General Linear Model Prediction
randomizationInference	Flexible Randomization-Based Inference
randomizeBE	Function to Create a Random List for Crossover Studies
randomizeR	Randomization for Clinical Trials
randomizr	Easy to Use Tools for Common Forms of Random Assignment
randomLCA	Random Effects Latent Class Analysis
randomNames	Function for Generating Random Names and a Dataset
random.polychor.pa	A Parallel Analysis With Polychoric Correlation Matrices
randomUniformForest	Random Uniform Forests for Classification, Regression and Unsupervised Learning
randstr	Generate Random Strings
randtests	Testing randomness in R
randtoolbox	Toolbox for Pseudo and Quasi Random Number Generation and RNG Tests
RandVar	Implementation of random variables
rangeBuilder	Occurrence Filtering and Generation of Species Range Polygons
rangeMapper	A Platform for the Study of Macroecology of Life History Traits
rangemodelR	Mid-Domain Effect and Species Richness Patterns
ranger	A Fast Implementation of Random Forests
RankAggreg	Weighted rank aggregation
Rankcluster	Model-Based Clustering for Multivariate Partial Ranking Data
rankdist	Distance Based Ranking Models
rankhazard	Rank-Hazard Plots
RankResponse	Ranking Responses in a Single Response Question or a Multiple Response Question
RANKS	Ranking of Nodes with Kernelized Score Functions
RANN	Fast Nearest Neighbour Search (Wraps Arya and Mount's ANN Library)
RANN.L1	Fast Nearest Neighbour Search (Wraps ANN Library) Using L1 Metric
RAP	Reversal Association Pattern
RapidPolygonLookup	Polygon lookup using kd trees
RAPIDR	Reliable Accurate Prenatal non-Invasive Diagnosis R package
RApiSerialize	R API Serialization
RAppArmor	Bindings to AppArmor and Security Related Linux Tools
rappdirs	Application Directories: Determine Where to Save Data, Caches, and Logs
rapport	A Report Templating System
rapportools	Miscellaneous (stats) helper functions with sane defaults for reporting
RArcInfo	Functions to import data from Arc/Info V7.x binary coverages
rareGE	Testing Gene-Environment Interaction for Rare Genetic Variants
rareNMtests	Ecological and biogeographical null model tests for comparing rarefaction curves
Rarity	Calculation of Rarity Indices for Species and Assemblages of Species
rARPACK	Solvers for Large Scale Eigenvalue and SVD Problems
RaschSampler	Rasch Sampler
rasclass	Supervised Raster Image Classification
rase	Range Ancestral State Estimation for Phylogeography and Comparative Analyses
raster	Geographic Data Analysis and Modeling
rasterVis	Visualization Methods for Raster Data
RateDistortion	Routines for Solving Rate-Distortion Problems
rateratio.test	Exact rate ratio test
raters	A Modification of Fleiss' Kappa in Case of Nominal and Ordinal Variables
rationalfun	Manipulation of Rational Functions
RAtmosphere	Standard Atmospheric profiles
rattle	Graphical User Interface for Data Mining in R
rAverage	Parameter estimation for the averaging model of Information Integration Theory
rAvis	Interface to the Bird-Watching Dataset Proyecto AVIS
rbamtools	Read and Write BAM (Binary Alignment) Files
rbefdata	BEFdata R package
rbenchmark	Benchmarking routine for R
RBerkeley	Oracle 'Berkeley DB' Interface for R
rBeta2009	The Beta Random Number and Dirichlet Random Vector Generating Functions
rbhl	Interface to the 'Biodiversity' 'Heritage' Library
RbioRXN	Process Rhea, KEGG, MetaCyc, Unipathway Biochemical Reaction Data
rbiouml	Interact with BioUML Server
rbison	Interface to the 'USGS' 'BISON' 'API'
Rbitcoin	R & bitcoin integration
rbitcoinchartsapi	R Package for the BitCoinCharts.com API
Rblpapi	R Interface to Bloomberg
rbmn	Handling Linear Gaussian Bayesian Networks
rbokeh	R Interface for Bokeh
Rborist	Extensible, Parallelizable Implementation of the Random Forest Algorithm
rbounds	Perform Rosenbaum bounds sensitivity tests for matched and unmatched data
RBPcurve	The Residual-Based Predictiveness Curve
rbugs	Fusing R and OpenBugs and Beyond
rbundler	Rbundler manages an application's dependencies systematically and repeatedly
rbvs	Ranking-Based Variable Selection
RCA	Relational Class Analysis
R.cache	Fast and Light-Weight Caching (Memoization) of Objects and Results to Speed Up Computations
RCALI	Calculation of the Integrated Flow of Particles between Polygons
rcanvec	Access and Plot CanVec and CanVec+ Data for Rapid Basemap Creation in Canada
Rcapture	Loglinear Models for Capture-Recapture Experiments
rCarto	This package builds maps with a full cartographic layout
RCassandra	R/Cassandra interface
rCBA	CBA Classifier for R
rcbalance	Large, Sparse Optimal Matching with Refined Covariate Balance
rcbsubset	Optimal Subset Matching with Refined Covariate Balance
rcdd	Computational Geometry
rcdk	rcdk - Interface to the CDK Libraries
rcdklibs	rcdklib - CDK libraries packaged for R
RCEIM	RCEIM - R Cross Entropy Inspired Method for Optimization
RcellData	Example Dataset for 'Rcell' Package
Rcereal	C++ Header Files of 'cereal'
Rcgmin	Conjugate Gradient Minimization of Nonlinear Functions
rchallenge	A Simple Data Science Challenge System
rchess	Chess Move, Generation/Validation, Piece Placement/ Movement, and Check/Checkmate/Stalemate Detection
Rchoice	Discrete Choice (Binary, Poisson and Ordered) Models with Random Parameters
rChoiceDialogs	rChoiceDialogs collection
rcicr	Reverse-Correlation Image-Classification Toolbox
RCircos	Circos 2D Track Plot
RClimMAWGEN	RClimMAWGEN (R Climate Index Multi-site Auto-regressive Weather GENeretor): a package to generate time series of climate indices from RMAWGEN generations
rClinicalCodes	R tools for integrating with the www.clinicalcodes.org repository
rclinicaltrials	Download Aggregate Trial Information and Results from ClinicalTrials.gov
RClone	Partially Clonal Populations Analysis
Rclusterpp	Linkable C++ clustering
rCMA	R-to-Java Interface for 'CMA-ES'
rcmdcheck	Run 'R CMD check' from 'R' and Capture Results
Rcmdr	R Commander
RcmdrMisc	R Commander Miscellaneous Functions
RcmdrPlugin.BCA	Rcmdr Plug-In for Business and Customer Analytics
RcmdrPlugin.coin	Rcmdr Coin Plug-In
RcmdrPlugin.depthTools	R commander Depth Tools Plug-In
RcmdrPlugin.DoE	R Commander Plugin for (industrial) Design of Experiments
RcmdrPlugin.doex	Rcmdr plugin for Stat 4309 course
RcmdrPlugin.EACSPIR	Plugin de R-Commander para el Manual 'EACSPIR'
RcmdrPlugin.EBM	Rcmdr Evidence Based Medicine Plug-in Package
RcmdrPlugin.EcoVirtual	Rcmdr EcoVirtual Plugin
RcmdrPlugin.epack	Rcmdr plugin for time series
RcmdrPlugin.Export	Export R Output to LaTeX or HTML
RcmdrPlugin.EZR	R Commander Plug-in for the EZR (Easy R) Package
RcmdrPlugin.FactoMineR	Graphical User Interface for FactoMineR
RcmdrPlugin.GWRM	R Commander Plug-in for Fitting Generalized Waring Regression Models
RcmdrPlugin.HH	Rcmdr Support for the HH package
RcmdrPlugin.IPSUR	An IPSUR Plugin for the R Commander
RcmdrPlugin.KMggplot2	An Rcmdr Plug-in for Kaplan-Meier Plots and Other Plots by Using the ggplot2 Package
RcmdrPlugin.lfstat	Rcmdr Plug-In for low flow analysis
RcmdrPlugin.MA	Graphical User Interface for Conducting Meta-Analyses in R
RcmdrPlugin.mosaic	Adds menu items to produce mosaic plots and assoc plots to Rcmdr
RcmdrPlugin.MPAStats	R Commander Plug-in for MPA Statistics
RcmdrPlugin.NMBU	R Commander Plug-in for University Level Applied Statistics
RcmdrPlugin.orloca	orloca Rcmdr Plug-in
RcmdrPlugin.plotByGroup	Rcmdr plots by group using lattice
RcmdrPlugin.pointG	Graphical POINT of view for questionnaire data Rcmdr Plug-In
RcmdrPlugin.qual	Rcmdr plugin for quality control course
RcmdrPlugin.RMTCJags	R MTC Jags Rcmdr Plugin
RcmdrPlugin.ROC	Rcmdr Receiver Operator Characteristic Plug-In PACKAGE
RcmdrPlugin.sampling	Tools for sampling in Official Statistical Surveys
RcmdrPlugin.SCDA	Rcmdr Plugin for Designing and Analyzing Single-case Experiments
RcmdrPlugin.seeg	Rcmdr Plugin for seeg
RcmdrPlugin.SLC	SLC Rcmdr Plug-in
RcmdrPlugin.SM	Rcmdr Sport Management Plug-In
RcmdrPlugin.sos	Efficiently search the R help pages
RcmdrPlugin.steepness	Steepness Rcmdr Plug-in
RcmdrPlugin.survival	R Commander Plug-in for the survival Package
RcmdrPlugin.TeachingDemos	Rcmdr Teaching Demos Plug-In
RcmdrPlugin.temis	Graphical Integrated Text Mining Solution
RcmdrPlugin.UCA	UCA Rcmdr Plug-in
RCMIP5	Tools for Manipulating and Summarizing CMIP5 Data
Rcolombos	Interface to Colombos Compendia using the Exposed REST API
RColorBrewer	ColorBrewer Palettes
RConics	Computations on Conics
rcorpora	A Collection of Small Text Corpora of Interesting Data
Rcplex	R interface to CPLEX
RCPmod	Regions of common profiles modelling with mixtures-of-experts
Rcpp	Seamless R and C++ Integration
Rcpp11	R and C++11
RcppAnnoy	'Rcpp' Bindings for 'Annoy', a Library for Approximate Nearest Neighbors
RcppAPT	Rcpp Interface to the APT Package Manager
RcppArmadillo	'Rcpp' Integration for the 'Armadillo' Templated Linear Algebra Library
RcppBDT	Rcpp bindings for the Boost Date_Time library
rcppbugs	R binding for cppbugs
RcppCCTZ	'Rcpp' Bindings for the 'CCTZ' Library
RcppClassic	Deprecated 'classic' Rcpp API
RcppClassicExamples	Examples using RcppClassic to interface R and C++
RcppCNPy	Read-Write Support for NumPy Files via Rcpp
RcppDE	Global Optimization by Differential Evolution in C++
RcppDL	Deep Learning Methods via Rcpp
RcppEigen	'Rcpp' Integration for the 'Eigen' Templated Linear Algebra Library
RcppExamples	Examples using 'Rcpp' to Interface R and C++
RcppFaddeeva	'Rcpp' Bindings for the 'Faddeeva' Package
RcppGSL	'Rcpp' Integration for 'GNU GSL' Vectors and Matrices
RcppOctave	Seamless Interface to Octave – And Matlab
RcppParallel	Parallel Programming Tools for 'Rcpp'
RcppProgress	An Interruptible Progress Bar with OpenMP Support for C++ in R Packages
RcppRedis	'Rcpp' Bindings for 'Redis' using the 'hiredis' Library
RcppRoll	Efficient Rolling / Windowed Operations
RcppShark	R Interface to the Shark Machine Learning Library
RcppSMC	Rcpp bindings for Sequential Monte Carlo
RcppStreams	Rcpp Integration of the Streamulus DSEL for Stream Processing
RcppTOML	'Rcpp' Bindings to Parser for Tom's Obvious Markup Language
RcppXts	Interface the xts API via Rcpp
RcppZiggurat	'Rcpp' Integration of Different "Ziggurat" Normal RNG Implementations
RCriteo	Loading Criteo Data into R
rcrossref	Client for Various 'CrossRef' 'APIs'
rcrypt	Symmetric File Encryption Using GPG
RCryptsy	Access to Cryptsy Crypto-Currency Exchange Public Information API via R
Rcsdp	R interface to the CSDP semidefinite programming library
rCUR	CUR decomposition package
RCurl	General Network (HTTP/FTP/...) Client Interface for R
Rd2roxygen	Convert Rd to Roxygen Documentation
rda	Shrunken Centroids Regularized Discriminant Analysis
RDataCanvas	Basic Runtime Support for Datacanvas.io
rdatacite	'DataCite' Client for 'OAI-PMH' Methods and their Search 'API'
rdatamarket	Data access API for DataMarket.com
rdd	Regression Discontinuity Estimation
rddtools	Toolbox for Regression Discontinuity Design ('RDD')
rDEA	Robust Data Envelopment Analysis (DEA) for R
rdetools	Relevant Dimension Estimation (RDE) in Feature Spaces
R.devices	Unified Handling of Graphics Devices
rdian	Client Library for The Guardian
RDIDQ	It perform Quality check on data
RDieHarder	R interface to the dieharder RNG test suite
Rdistance	Distance Sampling Analyses
RDML	Importing Real-Time Thermo Cycler (qPCR) Data from RDML Format Files
rDNA	R Bindings for the Discourse Network Analyzer
RDota	Data Analysis Toolbox for Dota2
Rdpack	Update and Manipulate Rd Documentation Objects
rdrobust	Robust Data-Driven Statistical Inference in Regression-Discontinuity Designs
rdrop2	Programmatic Interface to the 'Dropbox' API
rdryad	Access for Dryad Web Services
RDS	Respondent-Driven Sampling
Rdsdp	R Interface to DSDP Semidefinite Programming Library
Rdsm	Threads Environment for R
RDSTK	An R wrapper for the Data Science Toolkit API
rDVR	The rDVR package allows you to start stop and save a video server from within R
ReacTran	Reactive transport modelling in 1D, 2D and 3D
readbitmap	Simple Unified Interface to Read Bitmap Images (BMP,JPEG,PNG)
readBrukerFlexData	Reads Mass Spectrometry Data in Bruker *flex Format
readbulk	Read and Combine Multiple Data Files
reader	Suite of Functions to Flexibly Read Data from Files
readGenalex	Read, Write, Manipulate and Convert GenAlEx-Format Genotype Files
readMLData	Reading Machine Learning Benchmark Data Sets in Different Formats
readMzXmlData	Reads Mass Spectrometry Data in mzXML Format
readODS	Read ODS Files
readr	Read Tabular Data
readstata13	Import Stata Data Files
readxl	Read Excel Files
RealVAMS	Multivariate VAM Fitting
reams	Resampling-Based Adaptive Model Selection
Rearrangement	Monotonize Point and Interval Functional Estimates by Rearrangement
REBayes	Empirical Bayes Estimation and Inference in R
rebird	R Client for the eBird Database of Bird Observations
rebmix	Finite Mixture Modeling, Clustering & Classification
rebus	Build Regular Expressions in a Human Readable Way
rebus.base	Core Functionality for the 'rebus' Package
rebus.datetimes	Date and Time Extensions for the 'rebus' Package
rebus.numbers	Numeric Extensions for the 'rebus' Package
rebus.unicode	Unicode Extensions for the 'rebus' Package
RECA	Relevant Component Analysis for Supervised Distance Metric Learning
rechonest	R Interface to Echo Nest API
ReCiPa	Redundancy Control in Pathways databases
recluster	Ordination Methods for the Analysis of Beta-Diversity Indices
recoder	A Simple and Flexible Recoder
recommenderlab	Lab for Developing and Testing Recommender Algorithms
recommenderlabBX	Book-Crossing Dataset (BX) for 'recommenderlab'
recommenderlabJester	Jester Dataset for 'recommenderlab'
reconstructr	Session Reconstruction and Analysis
RecordLinkage	Record Linkage in R
Records	Record Values and Record Times
recosystem	Recommender System using Matrix Factorization
reda	Recurrent Event Data Analysis
REdaS	Companion Package to the Book 'R: Einführung durch angewandte Statistik'
redcapAPI	R Interface to REDCap
REDCapR	Interaction Between R and REDCap
RedditExtractoR	Reddit Data Extraction Toolkit
reddPrec	Reconstruction of Daily Data - Precipitation
redist	Markov Chain Monte Carlo Methods for Redistricting Simulation
redland	RDF Library Bindings in R
rEDM	Applications of Empirical Dynamic Modeling from Time Series
REEMtree	Regression Trees with Random Effects for Longitudinal (Panel) Data
ref	References for R
referenceIntervals	Reference Intervals
RefFreeEWAS	EWAS using Reference-Free DNA Methylation Mixture Deconvolution
refGenome	Gene and Splice Site Annotation Using Annotation Data from Ensembl and UCSC Genome Browsers
RefManageR	Straightforward 'BibTeX' and 'BibLaTeX' Bibliography Management
refset	Subsets with Reference Semantics
refund	Regression with Functional Data
refund.shiny	Interactive Plotting for Functional Data Analyses
refund.wave	Wavelet-Domain Regression with Functional Data
RegClust	Cluster analysis via regression coefficients
reGenotyper	Detecting Mislabeled Samples in Genetic Data
REGENT	Risk Estimation for Genetic and Environmental Traits
regexr	Readable Regular Expressions
registry	Infrastructure for R Package Registries
reglogit	Simulation-Based Regularized Logistic Regression
regpro	Nonparametric Regression
regress	Gaussian linear models with linear covariance structure
RegressionFactory	Expander Functions for Generating Full Gradient and Hessian from Single- and Multi-Slot Base Distributions
regRSM	Random Subspace Method (RSM) for Linear Regression
regsel	Variable Selection and Regression
regsubseq	Detect and Test Regular Sequences and Subsequences
regtest	Regression testing
rehh	Searching for Footprints of Selection using Haplotype Homozygosity Based Tests
rela	Item Analysis Package with Standard Errors
relaimpo	Relative importance of regressors in linear models
Relatedness	An Algorithm to Infer Relatedness
relations	Data Structures and Algorithms for Relations
relax	relax – R Editor for Literate Analysis and lateX
relaxnet	Relaxation of glmnet models (as in relaxed lasso, Meinshausen 2007)
relaxo	Relaxed Lasso
reldist	Relative Distribution Methods
relen	Compute Relative Entropy
relevent	Relational Event Models
Reliability	Functions for estimating parameters in software reliability models
ReliabilityTheory	Tools for Structural Reliability Analysis
reliaR	Package for some probability distributions
relimp	Relative Contribution of Effects in a Regression Model
relMix	Relationship Inference Based on Mixtures
relSim	Relative Simulator
relsurv	Relative Survival
RelValAnalysis	Relative Value Analysis
rem	Relational Event Models (REM)
remix	Remix your data
rEMM	Extensible Markov Model for Modelling Temporal Relationships Between Clusters
remMap	Regularized Multivariate Regression for Identifying Master Predictors
remote	Empirical Orthogonal Teleconnections in R
remoter	Remote R: Control a Remote R Session from a Local One
REndo	Fitting Linear Models with Endogenous Regressors when No External Instruments are Available
Renext	Renewal Method for Extreme Values Extrapolation
RenextGUI	GUI for Renext
rentrez	Entrez in R
Reol	R interface to the Encyclopedia of Life
ReorderCluster	Reordering the dendrogram according to the class labels
RepeatABEL	GWAS for Multiple Observations on Related Individuals
RepeatedHighDim	Global tests for expression data of high-dimensional sets of molecular features
repfdr	Replicability Analysis for Multiple Studies of High Dimension
repijson	Tools for Handling EpiJSON (Epidemiology Data) Files
replicatedpp2w	Two-Way ANOVA-Like Method to Analyze Replicated Point Patterns
replicationInterval	Replication Interval Functions
repmis	Miscellaneous Tools for Reproducible Research
repo	A Resource Manager for R Objects
repolr	Repeated Measures Proportional Odds Logistic Regression
ReporteRs	Microsoft Word, Microsoft PowerPoint and HTML Documents Generation
ReporteRsjars	External jars required for package ReporteRs
reportr	A General Message and Error Reporting System
reportRx	Tools for automatically generating reproducible clinical report
reports	Assist the Workflow of Writing Academic Articles and Other Reports
reporttools	Generate LaTeX Tables of Descriptive Statistics
REPPlab	R Interface to 'EPP-Lab', a Java Program for Exploratory Projection Pursuit
represent	Determine the representativity of two multidimensional data sets
reproducer	Reproduce Statistical Analyses and Meta-Analyses
REQS	R/EQS Interface
request	High Level 'HTTP' Client
rerddap	General Purpose Client for 'ERDDAP' Servers
reReg	Recurrent Event Regression
resample	Resampling Functions
resemble	Regression and Similarity Evaluation for Memory-Based Learning in Spectral Chemometrics
reservoir	Tools for Analysis, Design, and Operation of Water Supply Storages
reshape	Flexibly reshape data
reshape2	Flexibly Reshape Data: A Reboot of the Reshape Package
reshapeGUI	A GUI for the reshape2 and plyr packages
ResistorArray	electrical properties of resistor networks
ResourceSelection	Resource Selection (Probability) Functions for Use-Availability Data
RESS	Integrates R and Essentia
REST	RcmdrPlugin Easy Script Templates
restimizeapi	Functions for Working with the 'www.estimize.com' Web Services
restlos	Robust Estimation of Location and Scatter
restorepoint	Debugging with Restore Points
resumer	Build Resumes with R
rethinker	RethinkDB Client
retimes	Reaction Time Analysis
retistruct	Retinal Reconstruction Program
retrosheet	Import Professional Baseball Data from 'Retrosheet'
reutils	Talk to the NCBI EUtils
reval	Repeated Function Evaluation for Sensitivity Analysis
revealedPrefs	Revealed Preferences and Microeconomic Rationality
revealjs	R Markdown Format for 'reveal.js' Presentations
RevEcoR	Reverse Ecology Analysis on Microbiome
reweight	Adjustment of Survey Respondent Weights
rex	Friendly Regular Expressions
Rexperigen	R Interface to Experigen
rexpokit	R wrappers for EXPOKIT; other matrix functions
Rfacebook	Access to Facebook API via R
rFDSN	Get Seismic Data from the International Federation of Digital Seismograph Networks
rFerns	Random Ferns Classifier
RFGLS	Rapid Feasible Generalized Least Squares
RFgroove	Importance Measure and Selection for Groups of Variables with Random Forests
rfigshare	An R Interface to 'figshare'
R.filesets	Easy Handling of and Access to Files Organized in Structured Directories
RFinanceYJ	RFinanceYJ
rfishbase	R Interface to 'FishBase'
rfisheries	'Programmatic Interface to the 'openfisheries.org' API'
Rfit	Rank Estimation for Linear Models
RFLPtools	Tools to analyse RFLP data
RFmarkerDetector	Multivariate Analysis of Metabolomics Data using Random Forests
rfml	MarkLogic NoSQL Database Server in-Database Analytics for R
Rfmtool	Fuzzy Measure Tools for R
rfoaas	R Interface to FOAAS
RFOC	Graphics for Spherical Distributions and Earthquake Focal Mechanisms
RForcecom	Data Integration Feature for Force.com and Salesforce.com
rfordummies	Code Examples to Accompany the Book "R for Dummies"
rforensicbatwing	BATWING for calculating forensic trace-suspect match probabilities
RFormatter	R Source Code Formatter
rfPermute	Estimate Permutation p-Values for Random Forest Importance Metrics
RFreak	R/FrEAK interface
rfUtilities	Random Forests Model Selection and Performance Evaluation
RGA	A Google Analytics API Client
rgabriel	Gabriel Multiple Comparison Test and Plot the Confidence Interval on Barplot
rgam	Robust Generalized Additive Model
rGammaGamma	Gamma convolutions for methylation array background correction
rgbif	Interface to the Global 'Biodiversity' Information Facility 'API'
Rgbp	Hierarchical Modeling and Frequency Method Checking on Overdispersed Gaussian, Poisson, and Binomial Data
RGCCA	RGCCA and Sparse GCCA for multi-block data analysis
rgcvpack	R Interface for GCVPACK Fortran Package
rgdal	Bindings for the Geospatial Data Abstraction Library
RGENERATE	Tools To Generate Vector Time Series
RGENERATEPREC	Tools To Generate Daily-Precipitation Time Series
RGenetics	R packages for genetics research
rgenoud	R Version of GENetic Optimization Using Derivatives
rgeolocate	IP Address Geolocation
rgeos	Interface to Geometry Engine - Open Source (GEOS)
rgexf	Build, Import and Export GEXF Graph Files
rggobi	Interface between R and GGobi
rgho	Access WHO Global Health Observatory Data from R
RGIFT	Create quizzes in GIFT Format
rgl	3D Visualization Using OpenGL
rglobi	R Interface to Global Biotic Interactions
Rglpk	R/GNU Linear Programming Kit Interface
rglwidget	'rgl' in 'htmlwidgets' Framework
Rgnuplot	R Interface for Gnuplot
RGoogleAnalytics	R Wrapper for the Google Analytics API
RGoogleAnalyticsPremium	Unsampled Data in R for Google Analytics Premium Accounts
RgoogleMaps	Overlays on Google map tiles in R
rgp	R genetic programming framework
rgpui	UI for the RGP genetic programming framework
rgr	Applied Geochemistry EDA
RGraphics	Data and Functions from the Book R Graphics, Second Edition
rgrass7	Interface Between GRASS 7 Geographical Information System and R
rGroovy	Groovy Language Integration
RGtk2	R bindings for Gtk 2.8.0 and above
RGtk2Extras	Data frame editor and dialog making wrapper for RGtk2
RH2	DBI/RJDBC interface to h2 Database
rhandsontable	Interface to the 'Handsontable.js' Library
rHealthDataGov	Retrieve data sets from the HealthData.gov data API
rhosp	Side Effect Risks in Hospital : Simulation and Estimation
Rhpc	Permits *apply() Style Dispatch for 'HPC'
RhpcBLASctl	Control the Number of Threads on 'BLAS'
rHpcc	Interface between HPCC and R
RHRV	Heart rate variability analysis of ECG data
RHT	Regularized Hotelling's T-square Test for Pathway (Gene Set) Analysis
R.huge	Methods for Accessing Huge Amounts of Data [deprecated]
ri	ri: R package for performing randomization-based inference for experiments
RI2by2	Randomization inference for treatment effects on a binary outcome
riceware	A Diceware Passphrase Implementation
rich	Computes and compares species richnesses
RidgeFusion	R Package for Ridge Fusion in Statistical Learning
ridigbio	Interface to the iDigBio Data API
Ridit	Ridit Analysis (An extension of the Kruskal-Wallis Test.)
RIFS	Random Iterated Function System (RIFS)
RImageJROI	Read 'ImageJ' Region of Interest (ROI) Files
RImagePalette	Extract the Colors from Images
RImpala	Using Cloudera 'Impala' Through 'R'
rinat	Access iNaturalist data through APIs
rindex	Indexing for R
RInside	C++ Classes to Embed R in C++ Applications
RInSp	R Individual Specialization (RInSp)
rio	A Swiss-Army Knife for Data I/O
rioja	Analysis of Quaternary Science Data
Rip46	Utils for IP4 and IP6 Addresses
ripa	R Image Processing and Analysis
riskR	Risk Management
riskRegression	Risk Regression Models for Survival Analysis with Competing Risks
risksetROC	Riskset ROC curve estimation from censored survival data
riskSimul	Risk Quantification for Stock Portfolios under the T-Copula Model
RISmed	Download Content from NCBI Databases
Ritc	Isothermal Titration Calorimetry (ITC) Data Analysis
rite	The Right Editor to Write R
RItools	Randomization Inference Tools
riv	Robust instrumental variables estimator
rivernet	Read, Analyze and Plot River Networks
riverplot	Sankey or Ribbon Plots
rivervis	River Visualisation Tool
Rivivc	In vitro in vivo correlation linear level A
rivr	Steady and Unsteady Open-Channel Flow Computation
RJaCGH	Reversible Jump MCMC for the Analysis of CGH Arrays
rjade	A Clean, Whitespace-Sensitive Template Language for Writing HTML
RJafroc	Analysis of Data Acquired Using the Receiver Operating Characteristic Paradigm and Its Extensions
rjags	Bayesian Graphical Models using MCMC
rJava	Low-Level R to Java Interface
RJDBC	Provides access to databases through the JDBC interface
rje	Miscellaneous useful functions
r.jive	Perform JIVE Decompositions for Multi-Source Data
rJPSGCS	R-interface to Gene Drop Simulation from JPSGCS
Rjpstatdb	R interface of the Gateway to Advanced and User-friendly Statistics Service
RJSDMX	R Interface to SDMX Web Services
rjson	JSON for R
RJSONIO	Serialize R objects to JSON, JavaScript Object Notation
rjstat	Read and Write 'JSON-stat' Data Sets
rJython	R interface to Python via Jython
rkafka	Using Apache 'Kafka' Messaging Queue Through 'R'
rkafkajars	External Jars Required for Package 'rkafka'
RKEA	R/KEA Interface
RKEAjars	R/KEA Interface Jars
RKEEL	Using Keel in R Code
RKEELdata	Datasets from KEEL for it Use in RKEEL
RKEELjars	Java Executable .jar Files for 'RKEEL'
RKlout	Fetch Klout Scores for Twitter Users
rknn	Random KNN Classification and Regression
Rknots	Topological Analysis of Knotted Proteins, Biopolymers and 3D Structures
rkt	Mann-Kendall Test, Seasonal and Regional Kendall Tests
rkvo	Read Key/Value Pair Observations
Rlab	Functions and Datasets Required for ST370 class
Rlabkey	Data Exchange Between R and LabKey Server
rLakeAnalyzer	Package for the Analysis of Lake Physics
rleafmap	Interactive Maps with R and Leaflet
rlecuyer	R Interface to RNG with Multiple Streams
Rlibeemd	Ensemble Empirical Mode Decomposition (EEMD) and Its Complete Variant (CEEMDAN)
rLiDAR	LiDAR Data Processing and Visualization
rLindo	R Interface to LINDO API
Rlinkedin	Access to the LinkedIn API via R
rlist	A Toolbox for Non-Tabular Data Manipulation
rlm	Robust Fitting of Linear Model
Rlof	R Parallel Implementation of Local Outlier Factor(LOF)
RLogicalOps	Process Logical Operations
RLRsim	Exact (Restricted) Likelihood Ratio Tests for Mixed and Additive Models
rLTP	R Interface to LTP-Cloud Service
RLumModel	Modelling Ordinary Differential Equations Leading to Luminescence
RLumShiny	'Shiny' Applications for the R Package 'Luminescence'
RM2	Revenue Management and Pricing Package
rmaf	Refined Moving Average Filter
RMallow	Fit Multi-Modal Mallows' Models to ranking data
Rmalschains	Continuous Optimization using Memetic Algorithms with Local Search Chains (MA-LS-Chains) in R
RMark	R Code for Mark Analysis
rmarkdown	Dynamic Documents for R
rmatio	Read and Write Matlab Files
R.matlab	Read and Write MAT Files and Call MATLAB from Within R
RMAWGEN	Multi-site Auto-regressive Weather GENerator
RMC	Functions for fitting Markov models
rmcfs	The MCFS-ID Algorithm for Feature Selection and Interdependency Discovery
rmdformats	HTML Output Formats and Templates for 'rmarkdown' Documents
rmdshower	'R' 'Markdown' Format for 'shower' Presentations
RMediation	Mediation Analysis Confidence Intervals
rmeta	Meta-analysis
rmetasim	An Individual-Based Population Genetic Simulation Environment
R.methodsS3	S3 Methods Simplified
rmgarch	Multivariate GARCH Models
rminer	Data Mining Classification and Regression Methods
Rmisc	Rmisc: Ryan Miscellaneous
Rmixmod	An Interface for MIXMOD
RmixmodCombi	Combining Mixture Components for Clustering
RMixpanel	R API for Mixpanel
RMKdiscrete	Sundry Discrete Probability Distributions
rmngb	Miscellaneous Collection of Functions for Medical Data Analysis
RMOA	Connect R with MOA for Massive Online Analysis
RMOAjars	External jars required for package RMOA
RMongo	MongoDB Client for R
rmongodb	R-MongoDB driver
Rmonkey	A Survey Monkey R Client
Rmosek	The R-to-MOSEK Optimization Interface
rmp	Rounded Mixture Package. Performs Probability Mass Function Estimation with Nonparametric Mixtures of Rounded Kernels
Rmpfr	R MPFR - Multiple Precision Floating-Point Reliable
Rmpi	Interface (Wrapper) to MPI (Message-Passing Interface)
RMRAINGEN	RMRAINGEN (R Multi-site RAINfall GENeretor): a package to generate daily time series of rainfall from monthly mean values
rms	Regression Modeling Strategies
rms.gof	Root-mean-square goodness-of-fit test for simple null hypothesis
RMTstat	Distributions, Statistics and Tests derived from Random Matrix Theory
rmumps	Wrapper for MUMPS Library
RMySQL	Database Interface and 'MySQL' Driver for R
rnaseqWrapper	Wrapper for several R packages and scripts to automate RNA-seq analysis
RnavGraph	Using Graphs as a Navigational Infrastructure
RnavGraphImageData	Some image data used in the RnavGraph package demos
RNaviCell	Visualization of High-Throughput Data on Large-Scale Biological Networks
rnbn	Access NBN Data
RNCBIEUtilsLibs	EUtils libraries for use in the R environment
RNCEP	Obtain, Organize, and Visualize NCEP Weather Data
rncl	An Interface to the Nexus Class Library
RND	Risk Neutral Density Extraction Package
RndTexExams	Build Random Multiple Choice Exams
RNeo4j	Neo4j Driver for R
rneos	XML-RPC Interface to NEOS
rnetcarto	Fast Network Modularity and Roles Computation by Simulated Annealing (Rgraph C Library Wrapper for R)
RNetCDF	Interface to NetCDF Datasets
RNetLogo	Provides an Interface to the Agent-Based Modelling Platform NetLogo
RNewsflow	Tools for Analyzing Content Homogeneity and News Diffusion using Computational Text Analysis
RNeXML	Semantically Rich I/O for the 'NeXML' Format
rngSetSeed	Seeding the Default RNG with a Numeric Vector
rngtools	Utility functions for working with Random Number Generators
rngWELL	Toolbox for WELL Random Number Generators
rngwell19937	Random number generator WELL19937a with 53 or 32 bit output
RNiftyReg	Image Registration Using the NiftyReg Library
rNMF	Robust Nonnegative Matrix Factorization
rnn	Recurrent Neural Network
rnoaa	'NOAA' Weather Data from R
rNOMADS	An Interface to the NOAA Operational Model Archive and Distribution System
rnrfa	UK National River Flow Archive Data from R
ROAuth	R Interface For OAuth
RobAStBase	Robust Asymptotic Statistics
robCompositions	Robust Estimation for Compositional Data
robcor	Robust Correlations
robeth	R functions for robust statistics
robfilter	Robust Time Series Filters
RobLox	Optimally robust influence curves and estimators for location and scale
RobLoxBioC	Infinitesimally robust estimators for preprocessing omics data
robotstxt	A 'robots.txt' Parser and 'Webbot'/'Spider'/'Crawler' Permissions Checker
RobPer	Robust Periodogram and Periodicity Detection Methods
robreg3S	Three-Step Regression and Inference for Cellwise and Casewise Contamination
RobRex	Optimally robust influence curves for regression and scale
RobRSVD	Robust Regularized Singular Value Decomposition
RObsDat	Data Management for Hydrology and Beyond Using the Observations Data Model
robumeta	Robust Variance Meta-Regression
robust	Robust Library
RobustAFT	Truncated Maximum Likelihood Fit and Robust Accelerated Failure Time Regression for Gaussian and Log-Weibull Case
robustbase	Basic Robust Statistics
robustDA	Robust Mixture Discriminant Analysis
RobustEM	Robust Mixture Modeling Fitted via Spatial-EM Algorithm for Model-Based Clustering and Outlier Detection
robustfa	An Object Oriented Solution for Robust Factor Analysis
robustgam	Robust Estimation for Generalized Additive Models
robustHD	Robust Methods for High-Dimensional Data
robustlmm	Robust Linear Mixed Effects Models
robustloggamma	Robust estimation of the generalized log gamma model
RobustRankAggreg	Methods for robust rank aggregation
robustrao	An Extended Rao-Stirling Diversity Index to Handle Missing Data
robustreg	Robust Regression Functions
robustvarComp	Robust Estimation of Variance Component Models
robustX	eXperimental Functionality for Robust Statistics
ROC632	Construction of diagnostic or prognostic scoring system and internal validation of its discriminative capacities based on ROC curve and 0.633+ boostrap resampling
rocc	ROC based classification
rockchalk	Regression Estimation and Presentation
RockFab	Rock fabric and strain analysis tools
rococo	RObust rank COrrelation COefficient and test
ROCR	Visualizing the Performance of Scoring Classifiers
ROCS	Receiver Operating Characteristics Surface
ROCt	Time-Dependent ROC Curve Estimators and Expected Utility Functions
ROCwoGS	Non-parametric estimation of ROC curves without Gold Standard Test
RODBC	ODBC Database Access
RODBCDBI	Provides Access to Databases Through the ODBC Interface
RODBCext	Parameterized Queries Extension for RODBC
rodd	Optimal Discriminating Designs
RODM	R interface to Oracle Data Mining
ROI	R Optimization Infrastructure
ROI.plugin.glpk	ROI-plugin GLPK
ROI.plugin.quadprog	ROI-plugin quadprog
ROI.plugin.symphony	ROI-plugin symphony
rollply	Moving-Window Add-on for 'plyr'
ROMIplot	Plots Surfaces of Rates of Mortality Improvement
R.oo	R Object-Oriented Programming with or without References
Rook	Rook - a web server interface for R
RootsExtremaInflections	Finds roots, extrema and inflection points of a curve
rootSolve	Nonlinear Root Finding, Equilibrium and Steady-State Analysis of Ordinary Differential Equations
ropensecretsapi	R Package for the OpenSecrets.org API
ROpenWeatherMap	R Interface to OpenWeatherMap API
ROptEst	Optimally robust estimation
ROptEstOld	Optimally robust estimation - old version
ROptimizely	R Optimizely API
ROptRegTS	Optimally robust estimation for regression-type models
ror	Robust Ordinal Regression MCDA library
ROracle	OCI Based Oracle Database Interface for R
rorcid	Interface to the 'Orcid.org' 'API'
rorutadis	Robust Ordinal Regression UTADIS
ROSE	ROSE: Random Over-Sampling Examples
rosm	Plot Raster Map Tiles from Open Street Map and Other Sources
rotationForest	Fit and Deploy Rotation Forest Models
rotations	Tools for Working with Rotation Data
Rothermel	Rothermel fire spread model for R
rotl	Interface to the 'Open Tree of Life' API
roughrf	Roughened Random Forests for Binary Classification
RoughSetKnowledgeReduction	Simplification of Decision Tables using Rough Sets
RoughSets	Data Analysis Using Rough Set and Fuzzy Rough Set Theories
rowr	Row-Based Functions for R Objects
roxygen2	In-Source Documentation for R
royston	Royston's H Test: Multivariate Normality Test
RPANDA	Phylogenetic ANalyses of DiversificAtion
rpanel	Simple interactive controls for R using the tcltk library
rpart	Recursive Partitioning and Regression Trees
rpartitions	Code for integer partitioning
rpart.plot	Plot 'rpart' Models: An Enhanced Version of 'plot.rpart'
rpartScore	Classification trees for ordinal responses
rpart.utils	Tools for parsing and manipulating rpart objects, including generating machine readable rules
rpca	RobustPCA: Decompose a Matrix into Low-Rank and Sparse Components
rpcdsearch	Tools for the Construction of Clinical Code Lists for Primary Care Database Studies
RPCLR	RPCLR (Random-Penalized Conditional Logistic Regression)
Rpdb	Read, write, visualize and manipulate PDB files
rpdo	Pacific Decadal Oscillation Index
RPEnsemble	Random Projection Ensemble Classification
rpf	Response Probability Functions
rpg	Easy Interface to Advanced PostgreSQL Features
rphast	Interface to PHAST Software for Comparative Genomics
Rphylip	An R interface for PHYLIP
Rphylopars	Phylogenetic Comparative Tools for Missing Data and Within-Species Variation
rpivotTable	Build Powerful Pivot Tables and Dynamically Slice & Dice your Data
rPlant	Interface to the Agave API
rplexos	Read and Analyze 'PLEXOS' Solutions
rplos	Interface to the Search 'API' for 'PLoS' Journals
rplotengine	R as a plotting engine
RPMG	Graphical User Interface (GUI) for Interactive R Analysis Sessions
RPMM	Recursively Partitioned Mixture Model
rpnf	Point and Figure Package
Rpoppler	PDF Tools Based on Poppler
rportfolios	Random portfolio generation
RPostgreSQL	R interface to the PostgreSQL database system
rPowerSampleSize	Sample Size Computations Controlling the Type-II Generalized Family-Wise Error Rate
RPPairwiseDesign	Resolvable partially pairwise balanced design and Space-filling design via association scheme
RPPanalyzer	Reads, Annotates, and Normalizes Reverse Phase Protein Array Data
rPref	Database Preferences and Skyline Computation
RPresto	DBI Connector to Presto
rprime	Functions for Working with 'Eprime' Text Files
rprintf	Adaptive Builder for Formatted Strings
rprojroot	Finding Files in Project Subdirectories
RProtoBuf	R Interface to the Protocol Buffers API
rpsychi	Statistics for psychiatric research
rpubchem	rpubchem - Interface to the PubChem Collection
RPublica	ProPublica API Client
RPushbullet	R Interface to the Pushbullet Messaging Service
RPyGeo	ArcGIS Geoprocessing in R via Python
rPython	Package Allowing R to Call Python
RQDA	R-based Qualitative Data Analysis
rqPen	Penalized Quantile Regression
Rquake	Seismic Hypocenter Determination
RQuantLib	R Interface to the 'QuantLib' Library
rr	Statistical Methods for the Randomized Response Technique
Rramas	Matrix population models
rrBLUP	Ridge Regression and Other Kernels for Genomic Selection
rrBlupMethod6	Re-parametrization of RR-BLUP to allow for a fixed residual variance
rrcov	Scalable Robust Estimators with High Breakdown Point
rrcovHD	Robust Multivariate Methods for High Dimensional Data
rrcovNA	Scalable Robust Estimators with High Breakdown Point for Incomplete Data
Rrdrand	'DRNG' on Intel CPUs with the 'RdRand' Instruction for R
rredis	"Redis" Key/Value Database Client
rredlist	'IUCN' Red List Client
rrepast	Invoke 'Repast Simphony' Simulation Models
RRF	Regularized Random Forest
rriskDistributions	Fitting Distributions to Given Data or Known Quantiles
rrlda	Robust Regularized Linear Discriminant Analysis
RRNA	Secondary Structure Plotting for RNA
RRreg	Correlation and Regression Analyses for Randomized Response Data
R.rsp	Dynamic Generation of Scientific Reports
RRTCS	Randomized Response Techniques for Complex Surveys
RSA	Response Surface Analysis
RSADBE	Data related to the book "R Statistical Application Development by Example"
rsae	Robust Small Area Estimation
RSAGA	SAGA Geoprocessing and Terrain Analysis in R
Rsampletrees	Sampletrees Input/Output Processing
RSAP	SAP Netweaver RFC connector for R
rsatscan	Tools, Classes, and Methods for Interfacing with SaTScan Stand-Alone Software
rSCA	An R Package for Stepwise Cluster Analysis
rscala	Bi-Directional Interface Between R and Scala with Callbacks
rscimark	SciMark 2.0 Benchmark for Scientific and Numerical Computing
RSclient	Client for Rserve
rsconnect	Deployment Interface for R Markdown Documents and Shiny Applications
rscopus	Scopus Database API Interface
rscorecard	A Method to Download Department of Education College Scorecard Data
rscproxy	statconn: provides portable C-style interface to R (StatConnector)
RSDA	R to Symbolic Data Analysis
rsdepth	Ray Shooting Depth (i.e. RS Depth) functions for bivariate analysis
rsdmx	Tools for Reading SDMX Data and Metadata
RSeed	borenstein analysis
rseedcalc	Estimating the Proportion of Genetically Modified Seeds in Seedlots via Multinomial Group Testing
RSEIS	Seismic Time Series Analysis Tools
RSelenium	R bindings for Selenium WebDriver
rsem	Robust Structural Equation Modeling with Missing Data and Auxiliary Variables
Rserve	Binary R server
rSFA	Slow Feature Analysis in R
rsgcc	Gini methodology-based correlation and clustering analysis of microarray and RNA-Seq gene expression data
rsggm	Robust Sparse Gaussian Graphical Modeling via the Gamma-Divergence
RSGHB	Functions for Hierarchical Bayesian Estimation: A Flexible Approach
RSiena	Siena - Simulation Investigation for Empirical Network Analysis
rsig	Robust Signature Selection for Survival Outcomes
RsimMosaic	R Simple IMage Mosaic creation library
RSiteCatalyst	R Client for Adobe Analytics API V1.4
RSKC	Robust sparse K-means
rsm	Response-Surface Analysis
RSMET	Get Real-Time Meteorological Data in SMET Format
rsml	Plant Root System Markup Language (RSML) File Processing
RSNNS	Neural Networks in R using the Stuttgart Neural Network Simulator (SNNS)
rsnps	Get SNP (Single-Nucleotide Polymorphism) Data on the Web
RSNPset	Efficient Score Statistics for Genome-Wide SNP Set Analysis
RSocrata	Download or Upload 'Socrata' Data Sets
RSofia	Port of sofia-ml (http://code.google.com/p/sofia-ml/) to R
Rsolnp	General Non-Linear Optimization
Rsomoclu	Somoclu
rspa	Adapt Numerical Records to Fit (in)Equality Restrictions
rSPACE	Spatially-Explicit Power Analysis for Conservation and Ecology
RSpectra	Solvers for Large Scale Eigenvalue and SVD Problems
RSpincalc	Conversion Between Attitude Representations of DCM, Euler Angles, Quaternions, and Euler Vectors
RSPS	RNA-Seq Power Simulation
RSQLite	SQLite Interface for R
RSQLServer	SQL Server R Database Interface (DBI) and dplyr SQL Backend
Rssa	A Collection of Methods for Singular Spectrum Analysis
rstackdeque	Persistent Fast Amortized Stack and Queue Data Structures
rstan	R Interface to Stan
rstanarm	Bayesian Applied Regression Modeling via Stan
RStars	Access to the Digital Universe Data set API
RStata	A Bit of Glue Between R and Stata
rstatscn	R Interface for China National Data
rstiefel	Random orthonormal matrix generation on the Stiefel manifold
RStoolbox	Tools for Remote Sensing Data Analysis
RStorm	Simulate and Develop Streaming Processing in [R]
rstpm2	Generalized Survival Models
rstream	Streams of Random Numbers
rstudioapi	Safely Access the RStudio API
rsubgroup	Subgroup Discovery and Analytics
Rsundials	Suite of Nonlinear Differential Algebraic Equations Solvers in R
rsunlight	Interface to 'Sunlight' Foundation 'APIs'
Rsurrogate	Robust Estimation of the Proportion of Treatment Effect Explained by Surrogate Marker Information
RSurveillance	Design and Analysis of Disease Surveillance Activities
RSurvey	Analysis of Spatially Distributed Data
rsvd	Randomized Singular Value Decomposition
rsvg	Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays
RSvgDevice	An R SVG graphics device
RSVGTipsDevice	An R SVG graphics device with dynamic tips and hyperlinks
Rsymphony	SYMPHONY in R
rSymPy	R interface to SymPy computer algebra system
rtable	Tabular Reporting Functions
rTableICC	Random Generation of Contingency Tables
rtape	Manage and manipulate large collections of R objects stored as tape-like files
RTConnect	Tools for analyzing sales report files of iTunes Connect
RTDE	Robust Tail Dependence Estimation
rtdists	Response time distributions
rtematres	The rtematres API package
rTensor	Tools for Tensor Analysis and Decomposition
RTextTools	Automatic Text Classification via Supervised Learning
RTextureMetrics	Functions for calculation of texture metrics for Grey Level Co-occurrence Matrices
rtf	Rich Text Format (RTF) Output
rtfbs	Transcription Factor Binding Site Identification Tool
rticles	Article Formats for R Markdown
rtiff	Read and Write TIFF Files
rtimes	Client for New York Times 'APIs'
rtkore	STK++ Core Library Integration to R using Rcpp
rtkpp	STK++ Integration to R using Rcpp
RTOMO	Visualization for Seismic Tomography
rtop	Interpolation of Data with Variable Spatial Support
RTriangle	Triangle - A 2D Quality Mesh Generator and Delaunay Triangulator
rts	Raster Time Series Analysis
Rtsne	T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation
rtson	Typed JSON
Rttf2pt1	Package for ttf2pt1 program
Rtts	Convert Text into Speech
RtutoR	Tutorial App for Learning R
Rtwalk	The R Implementation of the 't-walk' MCMC Algorithm
rtype	A strong type system for R
Ruchardet	R package to detect character encoding
rucm	Implementation of Unobserved Components Model (UCM)
rugarch	Univariate GARCH Models
rUnemploymentData	Data and Functions for USA State and County Unemployment Data
RUnit	R Unit Test Framework
runittotestthat	Convert 'RUnit' Test Functions into 'testthat' Tests
Runiversal	Runiversal - Package for converting R objects to Java variables and XML
runjags	Interface Utilities, Model Templates, Parallel Computing Methods and Additional Distributions for MCMC Models in JAGS
Runuran	R Interface to the UNU.RAN Random Variate Generators
RunuranGUI	A GUI for the UNU.RAN random variate generators
rusda	Interface to USDA Databases
R.utils	Various Programming Utilities
ruv	Detect and Remove Unwanted Variation using Negative Controls
rv	Simulation-based random variable objects
RVAideMemoire	Diverse Basic Statistical and Graphical Functions
rvalues	R-Values for Ranking in High-Dimensional Settings
Rvcg	Manipulations of Triangular Meshes Based on the 'VCGLIB' API
rversions	Query 'R' Versions, Including 'r-release' and 'r-oldrel'
rvertnet	Search 'Vertnet', a 'Database' of Vertebrate Specimen Records
rvest	Easily Harvest (Scrape) Web Pages
RVFam	Rare Variants Association Analyses with Family Data
rvg	R Graphics Devices for Vector Graphics Output
rvgtest	Tools for Analyzing Non-Uniform Pseudo-Random Variate Generators
rvHPDT	Calling haplotype-based and variant-based pedigree disequilibrium test for rare variants in pedigrees
RVideoPoker	Play Video Poker with R
RViennaCL	ViennaCL C++ Header Files
Rvmmin	Variable Metric Nonlinear Function Minimization
RVowpalWabbit	R interface to the Vowpal Wabbit
RVPedigree	Methods for Family-Based Rare-Variant Genetic Association Tests
RVsharing	Probability of Sharing Rare Variants among Relatives
rvTDT	population control weighted rare-variants TDT
RVtests	Rare Variant Tests
Rwave	Time-Frequency Analysis of 1-D Signals
rWBclimate	A package for accessing World Bank climate data
RWBP	Detects spatial outliers using a Random Walk on Bipartite Graph
RWeather	R wrapper around the Yahoo! Weather, Google Weather and NOAA APIs
RWebLogo	plotting custom sequence logos
RWeka	R/Weka Interface
RWekajars	R/Weka Interface Jars
rwfec	R Wireless, Forward Error Correction
RWiener	Wiener process distribution functions
RWinEdt	R Interface to WinEdt
Rwinsteps	Running Winsteps in R
rwirelesscom	Basic Wireless Communications Simulation
rworldmap	Mapping Global Data
rworldxtra	Country boundaries at high resolution
rwt	Rice Wavelet Toolbox wrapper
rwunderground	R Interface to Weather Underground API
RxCEcolInf	R x C Ecological Inference With Optional Incorporation of Survey Information
RXKCD	Get XKCD comic from R
RXMCDA	Functions to Parse and Create XMCDA Files
RxnSim	Functions to Compute Chemical Reaction Similarity
RxODE	Facilities for Simulating from ODE-Based Models
RXshrink	Maximum Likelihood Shrinkage via Generalized Ridge or Least Angle Regression
Ryacas	R interface to the yacas computer algebra system
RYandexTranslate	R Interface to Yandex Translate API
RYoudaoTranslate	R package provide functions to translate English words into Chinese
ryouready	Companion to the Forthcoming Book - R you Ready?
rYoutheria	Access to the YouTheria mammal trait database
rysgran	Grain size analysis, textural classifications and distribution of unconsolidated sediments
Rz	GUI Tool for Data Management like SPSS or Stata
rzmq	R Bindings for ZeroMQ
0	
s20x	Functions for University of Auckland Course STATS 201/208 Data Analysis
s2dverification	Set of Common Tools for Forecast Verification
S2sls	Spatial Two Stage Least Squares Estimation
s4vd	Biclustering via Sparse Singular Value Decomposition Incorporating Stability Selection
Sabermetrics	Sabermetrics Functions for Baseball Analytics
sac	Semiparametric Analysis of Changepoint
saccades	Detection of Fixations in Eye-Tracking Data
SACCR	SA Counterparty Credit Risk under Basel III
SACOBRA	Self-Adjusting COBRA
sadists	Some Additional Distributions
sads	Maximum Likelihood Models for Species Abundance Distributions
sae	Small Area Estimation
sae2	Small Area Estimation: Time-series Models
saemix	Stochastic Approximation Expectation Maximization (SAEM) algorithm
SAENET	A Stacked Autoencoder Implementation with Interface to 'neuralnet'
saery	Small Area Estimation for Rao and Yu Model
saeSim	Simulation Tools for Small Area Estimation
SAFD	Statistical Analysis of Fuzzy Data
safeBinaryRegression	Safe Binary Regression
SafeQuant	A Toolbox for the Analysis of Proteomics Data
safi	Sensitivity Analysis for Functional Input
SAGA	Software for the Analysis of Genetic Architecture
SALES	Elastic Net and (Adaptive) Lasso Penalized Sparse Asymmetric Least Squares (SALES) and Coupled Sparse Asymmetric Least Squares (COSALES) using Coordinate Descent and Proximal Gradient Algorithms
SALTSampler	Efficient Sampling on the Simplex
SAM	Sparse Additive Modelling
SamplerCompare	A Framework for Comparing the Performance of MCMC Samplers
sampleSelection	Sample Selection Models
samplesize	Sample Size Calculation for Various t-Tests and Wilcoxon-Test
Sample.Size	Sample size calculation
samplesize4surveys	Sample Size Calculations for Complex Surveys
SampleSizeMeans	Sample size calculations for normal means
SampleSizeProportions	Calculating sample size requirements when estimating the difference between two binomial proportions
sampling	Survey Sampling
samplingbook	Survey Sampling Procedures
samplingEstimates	Sampling Estimates
SamplingStrata	Optimal Stratification of Sampling Frames for Multipurpose Sampling Surveys
samplingVarEst	Sampling Variance Estimation
sampSurf	Sampling Surface Simulation for Areal Sampling Methods
samr	SAM: Significance Analysis of Microarrays
SAMUR	Stochastic Augmentation of Matched Data Using Restriction Methods
SAMURAI	Sensitivity Analysis of a Meta-analysis with Unpublished but Registered Analytical Investigations
sand	Statistical Analysis of Network Data with R
sandwich	Robust Covariance Matrix Estimators
sanitizers	C/C++ source code to trigger Address and Undefined Behaviour Sanitizers
sankey	Sankey Diagrams
sanon	Stratified Analysis with Nonparametric Covariable Adjustment
sapa	Spectral Analysis for Physical Applications
SAPP	Statistical Analysis of Point Processes
sas7bdat	SAS Database Reader (experimental)
SAScii	Import ASCII files directly into R using only a SAS input script
SASmixed	Data sets from "SAS System for Mixed Models"
SASPECT	Significant AnalysiS of PEptide CounTs
SASxport	Read and Write 'SAS' 'XPORT' Files
satellite	Various Functions for Handling and Manipulating Remote Sensing Data
saturnin	Spanning Trees Used for Network Inference
SAVE	Bayesian Emulation, Calibration and Validation of Computer Models
saves	Fast load variables
saws	Small-Sample Adjustments for Wald tests Using Sandwich Estimators
sBF	Smooth Backfitting
sbgcop	Semiparametric Bayesian Gaussian copula estimation and imputation
sbioPN	sbioPN: Simulation of deterministic and stochastic spatial biochemical reaction networks using Petri Nets
sbmSDP	Semidefinite Programming for Fitting Block Models of Equal Block Sizes
SBRect	Detecting structural breaks using rectangle covering (non-parametric method)
SBSA	Simplified Bayesian Sensitivity Analysis
sca	Simple Component Analysis
scagnostics	Compute scagnostics - scatterplot diagnostics
Scale	Likert Type Questionnaire Item Analysis
scales	Scale Functions for Visualization
scalreg	Scaled sparse linear regression
scam	Shape Constrained Additive Models
scape	Statistical Catch-at-Age Plotting Environment
scar	Shape-Constrained Additive Regression: a Maximum Likelihood Approach
scaRabee	Optimization Toolkit for Pharmacokinetic-Pharmacodynamic Models
scatterD3	D3 JavaScript Scatterplot from R
scatterplot3d	3D Scatter Plot
SCBmeanfd	Simultaneous Confidence Bands for the Mean of Functional Data
scenario	Construct Reduced Trees with Predefined Nodal Structures
SCEPtER	Stellar CharactEristics Pisa Estimation gRid
SCEPtERbinary	Stellar CharactEristics Pisa Estimation gRid for Binary Systems
SCGLR	Supervised Component Generalized Linear Regression
SchemaOnRead	Automated Schema on Read
scholar	Analyse Citation Data from Google Scholar
schoolmath	Functions and datasets for math used in school
schoRsch	Tools for Analyzing Factorial Experiments
schumaker	Schumaker Shape-Preserving Spline
schwartz97	A package on the Schwartz two-factor commodity model
SCI	Standardized Climate Indices such as SPI, SRI or SPEI
scidb	An R Interface to SciDB
SciencesPo	A Tool Set for Analyzing Political Behavior Data
scio	Sparse Column-wise Inverse Operator
sciplot	Scientific Graphing Functions for Factorial Designs
SciViews	SciViews GUI API - Main package
sclero	Measure Growth Patterns and Align Sampling Spots in Photographs
SCMA	Single-Case Meta-Analysis
scmamp	Statistical Comparison of Multiple Algorithms in Multiple Problems
score	A Package to Score Behavioral Questionnaires
ScoreGGUM	Score Persons Using the Generalized Graded Unfolding Model
scorer	Quickly Score Models in Data Science and Machine Learning
SCORER2	SCORER 2.0: an algorithm for distinguishing parallel dimeric and trimeric coiled-coil sequences
scoring	Proper scoring rules
ScottKnott	The ScottKnott Clustering Algorithm
scout	Implements the Scout Method for Covariance-Regularized Regression
SCperf	Supply Chain Perform
ScrabbleScore	Calculates Scrabble score for strings
scrapeR	Tools for Scraping Data from HTML and XML Documents
ScreenClean	Screen and clean variable selection procedures
scrime	Analysis of High-Dimensional Categorical Data such as SNP Data
scriptests	Transcript-Based Unit Tests that are Easy to Create and Maintain
scrm	Simulating the Evolution of Biological Sequences
SCRSELECT	Performs Bayesian Variable Selection on the Covariates in a Semi-Competing Risks Model
SCRT	Single-Case Randomization Tests
scrubr	Clean Biological Occurrence Records
scrypt	scrypt key derivation functions for R
scs	Splitting Conic Solver
scuba	Diving Calculations and Decompression Models
SCVA	Single-Case Visual Analysis
sda	Shrinkage Discriminant Analysis and CAT Score Variable Selection
SDaA	Sampling: Design and Analysis
sdcMicro	Statistical Disclosure Control Methods for Anonymization of Microdata and Risk Estimation
sdcMicroGUI	Graphical User Interface for Package 'sdcMicro'
sdcTable	Methods for Statistical Disclosure Control in Tabular Data
sdcTarget	Statistical Disclosure Control Substitution Matrix Calculator
SDD	Serial Dependence Diagrams
SDDE	Shortcuts, Detours and Dead Ends (SDDE) Path Types in Genome Similarity Networks
sddpack	Semidiscrete Decomposition
sde	Simulation and Inference for Stochastic Differential Equations
sdef	Synthesizing List of Differentially Expressed Features
SDMTools	Species Distribution Modelling Tools: Tools for processing data associated with species distribution modelling exercises
sdmvspecies	Create Virtual Species for Species Distribution Modelling
sdnet	Soft Discretization-based Bayesian Network Inference
sdPrior	Scale-Dependent Hyperpriors in Structured Additive Distributional Regression
sdprisk	Measures of Risk for the Compound Poisson Risk Process with Diffusion
SDR	Subgroup Discovery Algorithms for R
sdtoolkit	Scenario Discovery Tools to Support Robust Decision Making
sdwd	Sparse Distance Weighted Discrimination
seacarb	Seawater Carbonate Chemistry
sealasso	Standard Error Adjusted Adaptive Lasso
searchable	Tools for Custom Searches / Subsets / Slices of Named R Objects
searchConsoleR	Google Search Console APIv3 R Client
SearchTrees	Spatial Search Trees
seas	Seasonal analysis and graphics, especially for climatology
SEAsic	Score Equity Assessment- summary index computation
season	Seasonal analysis of health data
seasonal	R Interface to X-13-ARIMA-SEATS
seawaveQ	U.S. Geological Survey seawaveQ model
SEchart	SEchart
SECP	Statistical Estimation of Cluster Parameters (SECP)
secr	Spatially Explicit Capture-Recapture
secrdesign	Sampling Design for Spatially Explicit Capture-Recapture
secrlinear	Spatially Explicit Capture-Recapture for Linear Habitats
seeclickfixr	Access Data from the SeeClickFix Web API
seedy	Simulation of Evolutionary and Epidemiological Dynamics
seeg	Statistics for Environmental Sciences, Engineering, and Geography
seem	Simulation of Ecological and Environmental Models
SEER2R	reading and writing SEER*STAT data files
SEERaBomb	SEER and Atomic Bomb Survivor Data Analysis Tools
seewave	Sound Analysis and Synthesis
seg	A set of tools for measuring spatial segregation
SegCorr	Detecting Correlated Genomic Regions
segmag	Determine Event Boundaries in Event Segmentation Experiments
segmented	Regression Models with Breakpoints/Changepoints Estimation
Segmentor3IsBack	A Fast Segmentation Algorithm
SEHmodel	Spatial Exposure-Hazard Model for Exposure and Impact Assessment on Exposed Individuals
seismic	Predict Information Cascade by Self-Exciting Point Process
seismicRoll	Fast Rolling Functions for Seismology using Rcpp
sejmRP	An Information About Deputies and Votings in Polish Diet from Seventh to Eighth Term of Office
Sejong	KoNLP static dictionaries and Sejong project resources
SEL	Semiparametric elicitation
selection	Correcting Biased Estimates Under Selection
selectiongain	A Tool for Calculation and Optimization of the Expected Gain from Multi-Stage Selection
selectiveInference	Tools for Post-Selection Inference
selectMeta	Estimation of Weight Functions in Meta Analysis
selectr	Translate CSS Selectors to XPath Expressions
selectspm	Select Point Pattern Models Based on Minimum Contrast, AIC and Goodness of Fit
SeleMix	Selective Editing via Mixture models
selfea	Select Features Reliably with Cohen's Effect Sizes
selfingTree	Genotype Probabilities in Intermediate Generations of Inbreeding Through Selfing
SelvarMix	Regularization for Variable Selection in Model-Based Clustering and Discriminant Analysis
sem	Structural Equation Models
semdiag	Structural equation modeling diagnostics
semGOF	Goodness-of-fit indexes for structural equation models
semiArtificial	Generator of Semi-Artificial Data
SemiCompRisks	Hierarchical Models for Parametric and Semi-Parametric Analyses of Semi-Competing Risks Data
SEMID	Identifiability of Linear Structural Equation Models
SemiMarkov	Multi-States Semi-Markov Models
SemiPar	Semiparametic Regression
SemiParBIVProbit	Semiparametric Copula Bivariate Regression Models
SemiParSampleSel	Semiparametric Sample Selection Modelling with Continuous or Discrete Response
semisupKernelPCA	Kernel PCA projection, and semi-supervised variant
SEMModComp	Model Comparisons for SEM
semPlot	Path diagrams and visual analysis of various SEM packages' output
semPLS	Structural Equation Modeling Using Partial Least Squares
semsfa	Semiparametric Estimation of Stochastic Frontier Models
semTools	Useful Tools for Structural Equation Modeling
sendmailR	send email using R
sendplot	Tool for sending interactive plots with tool-tip content
sensitivity	Sensitivity Analysis
sensitivity2x2xk	Sensitivity Analysis for 2x2xk Tables in Observational Studies
SensitivityCaseControl	Sensitivity Analysis for Case-Control Studies
sensitivitymv	Sensitivity Analysis in Observational Studies
sensitivitymw	Sensitivity analysis using weighted M-statistics
sensitivityPStrat	Principal Stratification Sensitivity Analysis Functions
SensMixed	Analysis of Sensory and Consumer Data in a Mixed Model Framework
SensoMineR	Sensory data analysis with R
sensory	Simultaneous Model-Based Clustering and Imputation via a Progressive Expectation-Maximization Algorithm
sensR	Thurstonian Models for Sensory Discrimination
SenSrivastava	Datasets from Sen & Srivastava
SensusR	Sensus Analytics
separationplot	Separation Plots
seqCBS	CN Profiling using Sequencing and CBS
seqDesign	Simulation and Group Sequential Monitoring of Randomized Two-Stage Treatment Efficacy Trials with Time-to-Event Endpoints
SeqFeatR	A Tool to Associate FASTA Sequences and Features
SeqGrapheR	Simple GUI for Graph Based Visualization of Cluster of DNA Sequence Reads
seqHMM	Hidden Markov Models for Life Sequences and Other Multivariate, Multichannel Categorical Time Series
seqinr	Biological Sequences Retrieval and Analysis
seqMeta	Meta-Analysis of Region-Based Tests of Rare DNA Variants
seqminer	Efficiently Read Sequence Data (VCF Format, BCF Format and METAL Format) into R
seqmon	Sequential Monitoring of Clinical Trials
seqPERM	Generates a permutation matrix based upon a sequence
seqRFLP	Simulation and visualization of restriction enzyme cutting pattern from DNA sequences
seqtest	Sequential Triangular Test
sequences	Generic and Biological Sequences
Sequential	Exact Sequential Analysis for Poisson and Binomial Data
sequenza	Copy Number Estimation from Tumor Genome Sequencing Data
serial	The Serial Interface Package
seriation	Infrastructure for Ordering Objects Using Seriation
seroincidence	Estimating Infection Rates from Serological Data
servr	A Simple HTTP Server to Serve Static Files or Dynamic Documents
sesem	Spatially explicit structural equation modeling
session	Functions for interacting with, saving and restoring R sessions
SetMethods	SetMethods: A Package Companion to "Set-Theoretic Methods for the Social Sciences"
SETPath	Spiked Eigenvalue Test for Pathway data
SetRank	Advanced Gene Set Enrichment Analysis
setRNG	Set (Normal) Random Number Generator and Seed
sets	Sets, Generalized Sets, Customizable Sets and Intervals
settings	Software Option Settings Manager for R
setwidth	Automatically Set the Width Option on Terminal Emulators
severity	Mayo's Post-data Severity Evaluation
sExtinct	Calculates the historic date of extinction given a series of sighting events
sfa	Stochastic Frontier Analysis
sfsmisc	Utilities from "Seminar fuer Statistik" ETH Zurich
sft	Functions for Systems Factorial Technology Analysis of Data
SGCS	Spatial Graph Based Clustering Summaries for Spatial Point Patterns
sgd	Stochastic Gradient Descent for Scalable Estimation
sgeostat	An Object-Oriented Framework for Geostatistical Modeling in S+
SGL	Fit a GLM (or cox model) with a combination of lasso and group lasso regularization
sglasso	Lasso Method for RCON(V,E) Models
sglOptim	Generic Sparse Group Lasso Solver
sglr	An R package for power and boundary calculations in pre-licensure vaccine trials using a sequential generalized likelihood ratio test
sgof	Multiple Hypothesis Testing
SGP	Student Growth Percentiles & Percentile Growth Trajectories
sGPCA	Sparse Generalized Principal Component Analysis
SGPdata	Exemplar Data Sets for SGP Analyses
sgPLS	Sparse Group Partial Least Square Methods
sgr	Sample Generation by Replacement
sgRSEA	Enrichment Analysis of CRISPR/Cas9 Knockout Screen Data
sgt	Skewed Generalized T Distribution Tree
shades	Simple Colour Manipulation
shape	Functions for plotting graphical shapes, colors
ShapeChange	Change-Point Estimation using Shape-Restricted Splines
shapefiles	Read and Write ESRI Shapefiles
shapeR	Collection and Analysis of Otolith Shape Data
shapes	Statistical Shape Analysis
ShapeSelectForest	Shape Selection for Landsat Time Series of Forest Dynamics
SharpeR	Statistical Significance of the Sharpe Ratio
sharpshootR	A Soil Survey Toolkit
sharx	Models and Data Sets for the Study of Species-Area Relationships
shazam	Immunoglobulin Somatic Hypermutation Analysis
SHELF	Tools to Support the Sheffield Elicitation Framework (SHELF)
shiny	Web Application Framework for R
shinyAce	Ace Editor Bindings for Shiny
shinybootstrap2	Bootstrap 2 Web Components for Use with Shiny
shinyBS	Twitter Bootstrap Components for Shiny
shinydashboard	Create Dashboards with 'Shiny'
shinyFiles	A Server-Side File System Viewer For Shiny
shinyjs	Perform Common JavaScript Operations in Shiny Apps using Plain R Code
shinyRGL	Shiny Wrappers for RGL
shinystan	Interactive Visual and Numerical Diagnostics and Posterior Analysis for Bayesian Models
shinythemes	Themes for Shiny
shinyTree	jsTree Bindings for Shiny
SHIP	SHrinkage covariance Incorporating Prior knowledge
shock	Slope Heuristic for Block-Diagonal Covariance Selection in High Dimensional Gaussian Graphical Models
shopifyr	An R Interface to the Shopify API
shotGroups	Analyze Shot Group Data
showtext	Using Fonts More Easily in R Graphs
showtextdb	Font Files for the 'showtext' Package
shp2graph	Convert a SpatialLinesDataFrame object to a "igraph-class" object
shrink	Global, Parameterwise and Joint Shrinkage Factor Estimation
Shrinkage	Several Shrinkage Effect-Size Estimators
ShrinkCovMat	Shrinkage Covariance Matrix Estimators
shuffle	The shuffle estimator for explainable variance
siar	Stable Isotope Analysis in R
SIBER	Stable Isotope Bayesian Ellipses in R
SID	Structural Intervention Distance
sideChannelAttack	Side Channel Attack
sidier	Substitution and Indel Distances to Infer Evolutionary Relationships
sievetest	Sieve test reporting functions
sig	Print Function Signatures
sigclust	Statistical Significance of Clustering
SightabilityModel	Wildlife Sightability Modeling
sigloc	Signal Location Estimation
signal	Signal Processing
signalHsmm	Predict Presence of Signal Peptides
signmedian.test	Perform Exact Sign Test and Asymptotic Sign Test in Large Samples
SigTree	Identify and Visualize Significantly Responsive Branches in a Phylogenetic Tree
SII	Calculate ANSI S3.5-1997 Speech Intelligibility Index
simba	A Collection of functions for similarity analysis of vegetation data
simboot	Simultaneous inference for diversity indices
simcausal	Simulating Longitudinal Data with Causal Inference Applications
SimComp	Simultaneous Comparisons for Multiple Endpoints
SimCorMultRes	Simulates Correlated Multinomial Responses
simctest	Safe Implementation of Monte Carlo Tests
SimDesign	Structure for Organizing Monte Carlo Simulation Designs
Sim.DiffProc	Simulation of Diffusion Processes
simecol	Simulation of Ecological (and Other) Dynamic Systems
simest	Constrained Single Index Model Estimation
simex	SIMEX- and MCSIMEX-Algorithm for measurement error models
simexaft	simexaft
simFrame	Simulation framework
SimHaz	Simulated Survival and Hazard Analysis for Time-Dependent Exposure
SimilarityMeasures	Trajectory Similarity Measures
Simile	Interact with Simile Models
SimInf	A Framework for Stochastic Disease Spread Simulations
simmer	Just Let it Simmer
simmr	A Stable Isotope Mixing Model
SIMMS	Subnetwork Integration for Multi-Modal Signatures
simMSM	Simulation of Event Histories for Multi-State Models
simone	Statistical Inference for MOdular NEtworks (SIMoNe)
simPH	Tools for Simulating and Plotting Quantities of Interest Estimated from Cox Proportional Hazards Models
simpleboot	Simple Bootstrap Routines
simplegraph	Simple Graph Data Types and Basic Algorithms
simpleNeural	An Easy to Use Multilayer Perceptron
simpleRCache	Simple R Cache
SimpleTable	Bayesian Inference and Sensitivity Analysis for Causal Effects from 2 x 2 and 2 x 2 x K Tables in the Presence of Unmeasured Confounding
simplexreg	Regression Analysis of Proportional Data Using Simplex Distribution
SimplicialCubature	Integration of Functions Over Simplices
simplr	Basic Symbolic Expression Simplification
simPop	Simulation of Synthetic Populations for Survey Data Considering Auxiliary Information
Simpsons	Detecting Simpson's Paradox
simr	Power Analysis for Generalised Linear Mixed Models by Simulation
SimRAD	Simulations to Predict the Number of RAD and GBS Loci
SimReg	Similarity Regression Functions
simrel	Linear Model Data Simulation and Design of Computer Experiments
simsalapar	Tools for Simulation Studies in Parallel
simsem	SIMulated Structural Equation Modeling
SimSeq	Nonparametric Simulation of RNA-Seq Data
simSummary	Simulation summary
simTool	Conduct Simulation Studies with a Minimal Amount of Source Code
SimuChemPC	Simulation process of 4 selection methods in predicting chemical potent compounds
SimultAnR	Correspondence and Simultaneous Analysis
SIN	A SINful Approach to Selection of Gaussian Graphical Markov Models
sinaplot	An Enhanced Chart for Simple and Truthful Representation of Single Observations over Multiple Classes
siplab	Spatial Individual-Plant Modelling
sirad	Functions for Calculating Daily Solar Radiation and Evapotranspiration
siRSM	Single-Index Response Surface Models
sirt	Supplementary Item Response Theory Models
SIS	Sure Independence Screening
sisal	Sequential Input Selection Algorithm
sisus	SISUS: Stable Isotope Sourcing using Sampling
sisVIVE	Some Invalid Some Valid Instrumental Variables Estimator
sitar	Super Imposition by Translation and Rotation Growth Curve Analysis
sitmo	Parallel Pseudo Random Number Generator (PPRNG) 'sitmo' Header Files
sitools	Format a number to a string with SI prefix
sivipm	Sensitivity Indices with Dependent Inputs
SixSigma	Six Sigma Tools for Quality Control and Improvement
SiZer	SiZer: Significant Zero Crossings
sjdbc	JDBC Driver Interface
sjmisc	Data Transformation and Labelled Data Utility Functions
sjPlot	Data Visualization for Statistics in Social Science
SKAT	SNP-Set (Sequence) Kernel Association Test
skatMeta	Efficient meta analysis for the SKAT test
skda	Sparse (Multicategory) Kernel Discriminant Analysis
skellam	Densities and Sampling for the Skellam Distribution
SkewHyperbolic	The Skew Hyperbolic Student t-Distribution
skewt	The Skewed Student-t Distribution
Skillings.Mack	The Skillings-Mack Test Statistic for Block Designs with Missing Observations
skmeans	Spherical k-Means Clustering
Sky	Canopy Openness Analyzer Package
sla	Two-Group Straight Line ANCOVA
slackr	Send messages, images, R objects and files to Slack.com channels/users
slam	Sparse Lightweight Arrays and Matrices
SLC	Slope and level change
sld	Estimation and Use of the Quantile-Based Skew Logistic Distribution
sleekts	4253H, Twice Smoothing
Sleuth2	Data Sets from Ramsey and Schafer's "Statistical Sleuth (2nd Ed)"
Sleuth3	Data Sets from Ramsey and Schafer's "Statistical Sleuth (3rd Ed)"
slfm	Tools for Fitting Sparse Latent Factor Model
SLHD	Maximin-Distance (Sliced) Latin Hypercube Designs
SLOPE	Sorted L1 Penalized Estimation (SLOPE)
slp	Discrete Prolate Spheroidal (Slepian) Sequence Regression Smoothers
sm	Smoothing methods for nonparametric regression and density estimation
smaa	Stochastic Multi-Criteria Acceptability Analysis
smac	Sparse Multi-category Angle-Based Large-Margin Classifiers
smacof	Multidimensional Scaling
smacpod	Statistical Methods for the Analysis of Case-Control Point Data
smallarea	Fits a Fay Herriot Model
smam	Statistical Modeling of Animal Movements
smart	Sparse Multivariate Analysis via Rank Transformation
SmarterPoland	Tools for Accessing Various Datasets Developed by the Foundation SmarterPoland.pl
smatr	(Standardised) Major Axis Estimation and Testing Routines
smbinning	Optimal Binning for Scoring Modeling
SMC	Sequential Monte Carlo (SMC) Algorithm
smcfcs	Multiple Imputation of Covariates by Substantive Model Compatible Fully Conditional Specification
smco	A simple Monte Carlo optimizer using adaptive coordinate sampling
SMCP	Smoothed minimax concave penalization (SMCP) method for genome-wide association studies
SMCRM	Data Sets for Statistical Methods in Customer Relationship Management by Kumar and Petersen (2012)
smcure	Fit Semiparametric Mixture Cure Models
smcUtils	Utility functions for sequential Monte Carlo
smdata	Data to accompany Smithson & Merkle, 2013
smdc	Document Similarity
smds	Symbolic Multidimensional Scaling
sme	Smoothing-splines Mixed-effects Models
smerc	Statistical Methods for Regional Counts
SMFI5	R functions and data from Chapter 5 of 'Statistical Methods for Financial Engineering'
smfsb	SMfSB 2e: Stochastic Modelling for Systems Biology, second edition
smint	Smooth Multivariate Interpolation for Gridded and Scattered Data
SMIR	Companion to Statistical Modelling in R
smirnov	Provides two taxonomic coefficients from E. S. Smirnov "Taxonomic analysis" (1969) book
SmithWilsonYieldCurve	Smith-Wilson Yield Curve Construction
SML	Statistical Machine Learning
SMNCensReg	Fitting Univariate Censored Regression Model Under the Family of Scale Mixture of Normal Distributions
smnet	Smoothing for Stream Network Data
smoof	Single and Multi-Objective Optimization Test Functions
smoother	Functions Relating to the Smoothing of Numerical Data
SmoothHazard	Fitting illness-death model for interval-censored data
smoothHR	Smooth Hazard Ratio Curves Taking a Reference Value
smoothie	Two-dimensional Field Smoothing
smoothmest	Smoothed M-estimators for 1-dimensional location
smoothSurv	Survival Regression with Smoothed Error Distribution
smoothtail	Smooth Estimation of GPD Shape Parameter
SMPracticals	Practicals for use with Davison (2003) Statistical Models
SMR	Externally Studentized Midrange Distribution
sms	Spatial Microsimulation
smss	Datasets for Agresti and Finlay's "Statistical Methods for the Social Sciences"
SMVar	Structural Model for variances
sn	The Skew-Normal and Skew-t Distributions
sna	Tools for Social Network Analysis
snapshot	Gadget N-body cosmological simulation code snapshot I/O utilities
SNFtool	Similarity Network Fusion
snht	Standard Normal Homogeneity Test
snipEM	Snipping methods for robust estimation and clustering
snn	Stabilized Nearest Neighbor Classifier
snow	Simple Network of Workstations
SnowballC	Snowball stemmers based on the C libstemmer UTF-8 library
snowboot	Network Analysis with Non-Parametric Methods that Emerge from Snowball and Bootstrap Sampling
snowfall	Easier cluster computing (based on snow)
snowFT	Fault Tolerant Simple Network of Workstations
snpar	Supplementary Non-parametric Statistics Methods
SNPassoc	SNPs-based whole genome association studies
snpEnrichment	SNPs Enrichment Analysis
snplist	Tools to Create Gene Sets
SNPmaxsel	Maximally selected statistics for SNP data
SNPMClust	Bivariate Gaussian Genotype Clustering and Calling for Illumina Microarrays
snp.plotter	snp.plotter
snpRF	Random Forest for SNPs to Prevent X-chromosome SNP Importance Bias
snpStatsWriter	Flexible writing of snpStats objects to flat files
SNPtools	Accessing, subsetting and plotting mouse SNPs
sns	Stochastic Newton Sampler (SNS)
SNscan	Scan Statistics in Social Networks
SNSequate	Standard and Nonstandard Statistical Models and Methods for Test Equating
SOAR	Memory management in R by delayed assignments
soc.ca	Specific Correspondence Analysis for the Social Sciences
SocialMediaLab	Tools for Collecting Social Media Data and Generating Networks for Analysis
SocialMediaMineR	A Social Media Search and Analytic Tool
SocialNetworks	Generates social networks based on distance
SocialPosition	Social Position Indicators Construction Toolbox
SOD	SOD for multidimensional scaling
SoDA	Functions and Examples for "Software for Data Analysis"
sodavis	SODA: Main and Interaction Effects Selection for Discriminant Analysis and Logistic Regression
SODC	Optimal Discriminant Clustering(ODC) and Sparse Optimal Discriminant Clustering(SODC)
sodium	A Modern and Easy-to-Use Crypto Library
Sofi	Interfaz interactiva con fines didacticos
softclassval	Soft Classification Performance Measures
SoftClustering	Soft Clustering Algorithms
softImpute	Matrix Completion via Iterative Soft-Thresholded SVD
soilDB	Soil Database Interface
soilphysics	Soil Physical Analysis
soilprofile	A package to consistently represent soil properties along a soil profile
SoilR	Models of Soil Organic Matter Decomposition
soil.spec	Soil Spectroscopy Tools and Reference Models
soiltexture	Functions for Soil Texture Plot, Classification and Transformation
soilwater	Implementation of Parametric Formulas for Soil Water Retention or Conductivity Curve
solaR	Radiation and Photovoltaic Systems
solarius	An R Interface to SOLAR
solarPos	Solar Position Algorithm for Solar Radiation Applications
solidearthtide	Solid Earth Tide Computation
SOLOMON	Parentage analysis
solr	General Purpose R Interface to Solr
solrium	General Purpose R Interface to 'Solr'
som	Self-Organizing Map
soma	General-Purpose Optimisation With the Self-Organising Migrating Algorithm
SOMbrero	SOM Bound to Realize Euclidean and Relational Outputs
somebm	some Brownian motions simulation functions
someKfwer	Controlling the Generalized Familywise Error Rate
someMTP	Some Multiple Testing Procedures
sommer	Solving Mixed Model Equations in R
somplot	Visualisation of hexagonal Kohonen maps
sonicLength	Estimating Abundance of Clones from DNA fragmentation data
soobench	Single Objective Optimization Benchmark Functions
SOPIE	Non-Parametric Estimation of the Off-Pulse Interval of a Pulsar
soql	Helps Make Socrata Open Data API Calls
SOR	Estimation using Sequential Offsetted Regression
SortableHTMLTables	Turns a data frame into an HTML file containing a sortable table
sortinghat	sortinghat
sorvi	Finnish Open Government Data Toolkit
sos	Search contributed R packages, sort by package
sos4R	An R client for the OGC Sensor Observation Service
sotkanet	Tools for Sotkanet Open Data Portal
soundecology	Soundscape Ecology
SoundexBR	Phonetic-Coding for Portuguese
SOUP	Stochastic Ordering Using Permutations (and Pairwise Comparisons)
source.gist	Read R code from a GitHub Gist
sourcetools	Tools for the Reading and Tokenization of R Code
SoyNAM	Soybean Nested Association Mapping Dataset
sp	Classes and Methods for Spatial Data
sp23design	Design and Simulation of seamless Phase II-III Clinical Trials
spa	Implements The Sequential Predictions Algorithm
SPA3G	SPA3G: R package for the method of Li and Cui (2012)
spaa	SPecies Association Analysis
space	Sparse PArtial Correlation Estimation
SPACECAP	A Program to Estimate Animal Abundance and Density using Bayesian Spatially-Explicit Capture-Recapture Models
spaceExt	Extension of SPACE
spacejam	Sparse conditional graph estimation with joint additive models
spacetime	Classes and Methods for Spatio-Temporal Data
spacodiR	Spatial and Phylogenetic Analysis of Community Diversity
spacom	Spatially Weighted Context Data for Multilevel Modelling
SpaDES	Develop and Run Spatially Explicit Discrete Event Simulation Models
spam	SPArse Matrix
spaMM	Mixed Models, Particularly Spatial GLMMs
spanel	Spatial Panel Data Models
spanr	Search Partition Analysis
SPAr	Perform rare variants association analysis based on summation of partition approaches
sparc	Semiparametric Generalized Linear Models
sparcl	Perform sparse hierarchical clustering and sparse k-means clustering
spareserver	Client Side Load Balancing
spark	Sparklines in the 'R' Terminal
sparkTable	Sparklines and Graphical Tables for TeX and HTML
sparktex	Generate LaTeX sparklines in R
SPARQL	SPARQL client
sparr	SPAtial Relative Risk
sparseBC	Sparse Biclustering of Transposable Data
sparsediscrim	Sparse and Regularized Discriminant Analysis
SparseFactorAnalysis	Scaling Count and Binary Data with Sparse Factor Analysis
SparseGrid	Sparse grid integration in R
sparseHessianFD	Numerical Estimation of Sparse Hessians
sparseLDA	Sparse Discriminant Analysis
SparseLearner	Sparse Learning Algorithms Using a LASSO-Type Penalty for Coefficient Estimation and Model Prediction
sparseLTSEigen	RcppEigen back end for sparse least trimmed squares regression
SparseM	Sparse Linear Algebra
sparseMVN	Multivariate Normal Functions for Sparse Covariance and Precision Matrices
sparsenet	Fit sparse linear regression models via nonconvex optimization
sparsereg	Sparse Bayesian Models for Regression, Subgroup Analysis, and Panel Data
sparseSEM	Sparse-aware Maximum Likelihood for Structural Equation Models
sparseSVM	Solution Paths of Sparse Linear Support Vector Machine with Lasso or ELastic-Net Regularization
SparseTSCGM	Sparse Time Series Chain Graphical Models
spartan	Simulation Parameter Analysis R Toolkit ApplicatioN: Spartan
spatcounts	Spatial count regression
spate	Spatio-Temporal Modeling of Large Data Using a Spectral SPDE Approach
spatgraphs	Graph Edge Computations for Spatial Point Patterns
spatial	Functions for Kriging and Point Pattern Analysis
spatialCovariance	Computation of Spatial Covariance Matrices for Data on Rectangles
spatialEco	Spatial Analysis and Modelling
SpatialEpi	Methods and Data for Spatial Epidemiology
SpatialExtremes	Modelling Spatial Extremes
spatialfil	Application of 2D Convolution Kernel Filters to Matrices or 3D Arrays
spatial.gev.bma	Hierarchical spatial generalized extreme value (GEV) modeling with Bayesian Model Averaging (BMA)
spatialkernel	Nonparameteric estimation of spatial segregation in a multivariate point process
spatialnbda	Performs spatial NBDA in a Bayesian context
SpatialNP	Multivariate nonparametric methods based on spatial signs and ranks
SpatialPack	Package for analysis of spatial data
SpatialPosition	Spatial Position Models
spatialprobit	Spatial Probit Models
spatialsegregation	Segregation measures for multitype spatial point patterns
spatialTailDep	Estimation of spatial tail dependence models
spatial.tools	R functions for working with spatial data
SpatialTools	Tools for Spatial Data Analysis
SpatialVx	Spatial Forecast Verification
SpatioTemporal	Spatio-Temporal Model Estimation
SpatPCA	Regularized Principal Component Analysis for Spatial Data
spatstat	Spatial Point Pattern Analysis, Model-Fitting, Simulation, Tests
spatsurv	Bayesian Spatial Survival Analysis with Parametric Proportional Hazards Models
spBayes	Univariate and Multivariate Spatial-temporal Modeling
spBayesSurv	Bayesian Modeling and Analysis of Spatially Correlated Survival Data
spc	Statistical Process Control – Collection of Some Useful Functions
spcadjust	Functions for Calibrating Control Charts
SPCALDA	A New Reduced-Rank Linear Discriminant Analysis Method
spcosa	Spatial Coverage Sampling and Random Sampling from Compact Geographical Strata
spcov	Sparse Estimation of a Covariance Matrix
spcr	Sparse principal component regression
spd	Semi Parametric Distribution
spdep	Spatial Dependence: Weighting Schemes, Statistics and Models
spduration	Split-Population Duration (Cure) Regression
spdynmod	Spatio-Dynamic Wetland Plant Communities Model
spe	Stochastic Proximity Embedding
speaq	Tools for Nuclear Magnetic Resonance Spectrum Alignment and Quantitative Analysis
speccalt	Alternative spectral clustering, with automatic estimation of k
SpecHelpers	Spectroscopy Related Utilities
SPECIES	Statistical package for species richness estimation
speciesgeocodeR	Prepare Species Distributions for the Use in Phylogenetic Analyses
SpeciesMix	Fit Mixtures of Archetype species
specificity	Specificity of personality trait-outcome (or trait-trait) associations
specmine	Metabolomics and Spectral Data Analysis and Mining
SpecsVerification	Forecast Verification Routines for the SPECS FP7 Project
spectralGP	Approximate Gaussian Processes Using the Fourier Basis
spectral.methods	Singular Spectrum Analysis (SSA) Tools for Time Series Analysis
spectrino	Spectra Visualization, Organizer and Data Preparation
speedglm	Fitting Linear and Generalized Linear Models to Large Data Sets
speff2trial	Semiparametric efficient estimation for a two-sample treatment effect
SPEI	Calculation of the Standardised Precipitation-Evapotranspiration Index
sperich	Auxiliary Functions to Estimate Centers of Biodiversity
sperrorest	Spatial Error Estimation and Variable Importance
spfrontier	Spatial Stochastic Frontier models estimation
spgrass6	Interface Between GRASS 6+ Geographical Information System and R
spgs	Statistical Patterns in Genomic Sequences
spgwr	Geographically Weighted Regression
sphereplot	Spherical plotting
SphericalCubature	Numerical integration over spheres and balls in n-dimensions; multivariate polar coordinates
SphericalK	Spherical K-Function
SpherWave	Spherical Wavelets and SW-based Spatially Adaptive Methods
sphet	Estimation of spatial autoregressive models with and without heteroskedastic innovations
spi	Compute SPI index
SPIAssay	A genetic-based assay for the identification of cell lines
spider	Species Identity and Evolution in R
spiders	Fits Predator Preferences Model
spikeslab	Prediction and variable selection using spike and slab regression
spikeSlabGAM	Bayesian Variable Selection and Model Choice for Generalized Additive Mixed Models
SPIn	Simulation-efficient Shortest Probability Intervals
spinyReg	Sparse Generative Model and Its EM Algorithm
splancs	Spatial and Space-Time Point Pattern Analysis
splitstackshape	Stack and Reshape Datasets After Splitting Concatenated Values
splm	Econometric Models for Spatial Panel Data
spls	Sparse Partial Least Squares (SPLS) Regression and Classification
splus2R	Supplemental S-PLUS functionality in R
splusTimeDate	Times and Dates from S-PLUS
splusTimeSeries	Time Series from S-PLUS
spm12r	Wrapper Functions for SPM (Statistical Parametric Mapping) Version 12 from the Wellcome Trust Centre for Neuroimaging
spMC	Continuous-Lag Spatial Markov Chains
SPmlficmcm	Semiparametric Maximum Likelihood Method for Interactions Gene-Environment in Case-Mother Control-Mother Designs
spnet	Plotting (Social) Networks on Maps
spocc	Interface to Species Occurrence Data Sources
spoccutils	Utilities for Use with 'spocc'
SPODT	Spatial Oblique Decision Tree
sporm	Semiparametric proportional odds rate model
SportsAnalytics	Infrastructure for Sports Analytics
SPOT	Sequential Parameter Optimization Toolbox
spray	Sparse Arrays and Multivariate Polynomials
SPREDA	Statistical Package for Reliability Data Analysis
sprex	Calculate Species Richness and Extrapolation Metrics
sprint	Simple Parallel R INTerface
sprinter	Framework for Screening Prognostic Interactions
sprintfr	An Easy Interface to String Formatting
sprm	Sparse and Non-Sparse Partial Robust M Regression and Classification
sprsmdl	Sparse modeling toolkit
SPRT	Wald's Sequential Probability Ratio Test
spsann	Optimization of Sample Configurations using Spatial Simulated Annealing
spsi	Shape-Preserving Uni-Variate and Bi-Variate Spline Interpolation
SPSL	Site Percolation on Square Lattice (SPSL)
spsmooth	spsmooth: An Extension Package for 'mgcv'
spsurvey	Spatial Survey Design and Analysis
spt	Sierpinski Pedal Triangle
spTDyn	Spatially Varying and Spatio-Temporal Dynamic Linear Models
spTest	Nonparametric Hypothesis Tests of Isotropy and Symmetry
spThin	Functions for Spatial Thinning of Species Occurrence Records for Use in Ecological Models
spTimer	Spatio-Temporal Bayesian Modelling Using R
sptm	SemiParametric Transformation Model Methods
spuRs	Functions and Datasets for "Introduction to Scientific Programming and Simulation Using R"
SQDA	Sparse Quadratic Discriminant Analysis
sqldf	Perform SQL Selects on R Data Frames
sqliter	Connection wrapper to SQLite databases
sqlutils	Utilities for working with SQL files
SQN	subset quantile normalization
SQUAREM	Squared extrapolation methods for accelerating fixed-point iterations
squash	Color-Based Plots for Multivariate Visualization
sra	Selection Response Analysis
SRCS	Statistical Ranking Color Scheme for Multiple Pairwise Comparisons
srd	Draws Scaled Rectangle Diagrams
sROC	Nonparametric Smooth ROC Curves for Continuous Data
SRRS	The Stepwise Response Refinement Screener (SRRS)
srvyr	'dplyr'-Like Syntax for Summary Statistics of Survey Data
ss3sim	Fisheries Stock Assessment Simulation Testing with Stock Synthesis
ssa	Simultaneous Signal Analysis
ssanv	Sample Size Adjusted for Nonadherence or Variability of Input Parameters
sscor	Robust Correlation Estimation and Testing Based on Spatial Signs
ssd	Sample Size Determination (SSD) for Unordered Categorical Data
SSDforR	Functions to Analyze Single System Data
SSDM	Stacked Species Distribution Modelling
sSDR	Tools Developed for Structured Sufficient Dimension Reduction (sSDR)
ssfa	Spatial Stochastic Frontier Analysis
ssfit	Fitting of parametric models using summary statistics
ssh.utils	Local and remote system commands with output and error capture
ssize.fdr	Sample Size Calculations for Microarray Experiments
ssizeRNA	Sample Size Calculation for RNA-Seq Experimental Design
ssmrob	Robust estimation and inference in sample selection models
SSN	Spatial Modeling on Stream Networks
sspline	Smoothing Splines on the Sphere
sspse	Estimating Hidden Population Size using Respondent Driven Sampling Data
SSrat	Two-dimensional sociometric status determination with rating scales
SSRMST	Sample Size Calculation using Restricted Mean Survival Time
sss	Tools for importing files in the triple-s .(Standard Survey Structure) format
SSsimple	State space models
ssvd	Sparse SVD
ssym	Fitting Semi-Parametric log-Symmetric Regression Models
st	Shrinkage t Statistic and Correlation-Adjusted t-Score
stabledist	Stable Distribution Functions
StableEstim	Estimate the 4 parameters of stable law using different methods
stabs	Stability Selection with Error Control
Stack	Stylized concatenation of data.frames or ffdfs
stackoverflow	Stack Overflow's Greatest Hits
stacomirtools	stacomi ODBC connection class
stagePop	Modelling the Population Dynamics of a Stage-Structured Species in Continuous Time
stam	Spatio-Temporal Analysis and Modelling
StAMPP	Statistical Analysis of Mixed Ploidy Populations
STAND	Statistical Analysis of Non-Detects
StandardizeText	Standardize Text
StanHeaders	C++ Header Files for Stan
STAR	Spike Train Analysis with R
stargazer	Well-Formatted Regression and Summary Statistics Tables
starma	Modelling Space Time AutoRegressive Moving Average (STARMA) Processes
startupmsg	Utilities for start-up messages
stashR	A Set of Tools for Administering SHared Repositories
Stat2Data	Datasets for Stat2
statar	Tools Inspired by Stata to Manipulate Tabular Data
statcheck	Extract Statistics from Articles and Recompute P Values
StatDA	Statistical Analysis for Environmental Data
StatDataML	Read and Write StatDataML Files
statebins	U.S. State Cartogram Heatmaps in R; an Alternative to Choropleth Maps for USA States
statfi	statfi R tools
stationaRy	Get Hourly Meteorological Data from Global Stations
StatMatch	Statistical Matching
StatMeasures	Easy Data Manipulation, Data Quality and Statistical Checks
StatMethRank	Statistical Methods for Ranking Data
statmod	Statistical Modeling
statnet	Software Tools for the Statistical Analysis of Network Data
statnet.common	Common R Scripts and Utilities Used by the Statnet Project Software
statnetWeb	A Graphical User Interface for Network Modeling with 'Statnet'
Statomica	Statomica utility package
staTools	Statistical Tools for Social Network Analysis
StatRank	Statistical Rank Aggregation: Inference, Evaluation, and Visualization
stdReg	Regression Standardization
steadyICA	ICA and Tests of Independence via Multivariate Distance Covariance
steepness	Testing Steepness of Dominance Hierarchies
SteinIV	Semi-Parametric Stein-Like Estimator with Instrumental Variables
stellaR	stellar evolution tracks and isochrones
Stem	Spatio-temporal models in R
STEPCAM	ABC-SMC Inference of STEPCAM
stepp	Subpopulation Treatment Effect Pattern Plot (STEPP)
stepPlr	L2 penalized logistic regression with a stepwise variable selection
stepR	Fitting Step-Functions
stepwise	Stepwise detection of recombination breakpoints
StereoMorph	Stereo Camera Calibration and Reconstruction
stheoreme	Klimontovich's S-Theorem Algorithm Implementation and Data Preparation Tools
STI	Calculation of the Standardized Temperature Index
stilt	Separable Gaussian Process Interpolation (Emulation)
stima	Simultaneous Threshold Interaction Modeling Algorithm
stinepack	Stineman, a consistently well behaved method of interpolation
stlplus	Enhanced Seasonal Decomposition of Time Series by Loess
stm	Estimation of the Structural Topic Model
stmBrowser	Structural Topic Model Browser
stmCorrViz	A Tool for Structural Topic Model Visualizations
STMedianPolish	Spatio-Temporal Median Polish
StMoMo	Stochastic Mortality Modelling
StMoSim	Plots a QQ-Norm Plot with several Gaussian simulations
stocc	Fit a Spatial Occupancy Model via Gibbs Sampling
stochprofML	Stochastic Profiling using Maximum Likelihood Estimation
stochvol	Efficient Bayesian Inference for Stochastic Volatility (SV) Models
StockChina	Real-Time Stock Price & Volume in China Market
stockPortfolio	Build stock models and analyze stock portfolios
stocks	Fast Functions for Stock Market Analysis
stoichcalc	R Functions for Solving Stoichiometric Equations
Storm	Write Storm Bolts in R using the Storm Multi-Language Protocol
storr	Simple Key Value Stores
stosim	Stochastic Simulator for Reliability Modeling of Repairable Systems
STPGA	Selection of Training Populations by Genetic Algorithm
stplanr	Sustainable Transport Planning
stpm	Stochastic Model for Analysis of Longitudinal Data
stpp	Space-Time Point Pattern simulation, visualisation and analysis
stppResid	Perform residual analysis on space-time point process models
StrainRanking	Ranking of pathogen strains
strap	Stratigraphic Tree Analysis for Palaeontology
strataG	Summaries and Population Structure Analyses of Haplotypic and Genotypic Data
stratification	Univariate Stratification of Survey Populations
stratigraph	Toolkit for the plotting and analysis of stratigraphic and palaeontological data
StratSel	Strategic Selection Estimator
straweib	Stratified Weibull Regression Model
stream	Infrastructure for Data Stream Mining
StreamMetabolism	Calculate Single Station Metabolism from Diurnal Oxygen Curves
streamMOA	Interface for MOA Stream Clustering Algorithms
streamR	Access to Twitter Streaming API via R
stremo	Functions to help the process of learning structural equation modelling
stressr	Fetch and plot financial stress index and component data
StressStrength	Computation and estimation of reliability of stress-strength models
stringdist	Approximate String Matching and String Distance Functions
stringgaussnet	PPI and Gaussian Network Construction from Transcriptomic Analysis Results Integrating a Multilevel Factor
stringi	Character String Processing Facilities
stringr	Simple, Consistent Wrappers for Common String Operations
stripless	Structured Trellis Displays Without Strips for Lattice Graphics
strptimer	An Easy Interface to Time Formatting
strucchange	Testing, Monitoring, and Dating Structural Changes
structSSI	Multiple Testing for Hypotheses with Hierarchical or Group Structure
strum	STRUctural Modeling of Latent Variables for General Pedigree
strvalidator	Process Control and Internal Validation of Forensic STR Kits
stsm	Structural Time Series Models
stsm.class	Class and Methods for Structural Time Series Models
stubthat	Stubbing Framework for R
stylo	Functions for a Variety of Stylometric Analyses
SubCultCon	Maximum-Likelihood Cultural Consensus Analysis with Sub-Cultures
subgroup	Methods for exploring treatment effect heterogeneity in subgroup analysis of clinical trials
SubLasso	Gene selection using Lasso for Microarray data with user-defined genes fixed in model
SubpathwayGMir	Identify Metabolic Subpathways Mediated by MicroRNAs
SubpathwayLNCE	Identify Signal Subpathways Competitively Regulated by LncRNAs Based on ceRNA Theory
subplex	Unconstrained Optimization using the Subplex Algorithm
subrank	Computes Copula using Ranks and Subsampling
subscore	SubScore Computing Functions in Classical Test Theory
subselect	Selecting Variable Subsets
subsemble	An Ensemble Method for Combining Subset-Specific Algorithm Fits
subspace	Interface to OpenSubspace
subtype	Cluster analysis to find molecular subtypes and their assessment
sudoku	Sudoku Puzzle Generator and Solver
sudokuAlt	Tools for Making and Spoiling Sudoku Games
SUE	Subsampling method
summarytools	Dataframe Summaries, Frequency Tables and Numerical Summaries with Customizable Output
Sunder	Quantification of the effect of geographic versus environmental isolation on genetic differentiation
SunterSampling	Sunter's sampling design
supclust	Supervised Clustering of Predictor Variables such as Genes
supcluster	Supervised Cluster Analysis
superbiclust	Generating Robust Biclusters from a Bicluster Set (Ensemble Biclustering)
superdiag	R Code for Testing Markov Chain Nonconvergence
SuperExactTest	Exact Test and Visualization of Multi-Set Intersections
SuperLearner	Super Learner Prediction
superMDS	Implements the supervised multidimensional scaling (superMDS) proposal of Witten and Tibshirani (2011)
superpc	Supervised principal components
SuppDists	Supplementary Distributions
support.BWS	Basic Functions for Supporting an Implementation of Best-Worst Scaling
support.CEs	Basic Functions for Supporting an Implementation of Choice Experiments
surface	Fitting Hansen Models to Investigate Convergent Evolution
Surrogate	Evaluation of Surrogate Endpoints in Clinical Trials
suRtex	LaTeX descriptive statistic reporting for survey data
surv2sampleComp	Inference for model-free between-group parameters for censored survival data
survAccuracyMeasures	Estimate accuracy measures for risk prediction markers from survival data
survAUC	Estimators of prediction accuracy for time-to-event data
survC1	C-statistics for risk prediction models with censored survival data
SurvCorr	Correlation of Bivariate Survival Times
surveillance	Temporal and Spatio-Temporal Modeling and Monitoring of Epidemic Phenomena
survexp.fr	Relative survival, AER and SMR based on French death rates
survey	analysis of complex survey samples
surveydata	Tools to manipulate survey data
surveyeditor	Generate a Survey that can be Completed by Survey Respondents
surveyoutliers	Helps Manage Outliers in Sample Surveys
surveyplanning	Survey Planning Tools
Survgini	The Gini concentration test for survival data
survIDINRI	IDI and NRI for comparing competing risk prediction models with censored survival data
survival	Survival Analysis
survivalMPL	Penalised Maximum Likelihood for Survival Analysis Models
survivalROC	Time-dependent ROC curve estimation from censored survival data
survJamda	Survival Prediction by Joint Analysis of Microarray Gene Expression Data
survJamda.data	Data for Package 'survJambda'
SurvLong	Analysis of Proportional Hazards Model with Sparse Longitudinal Covariates
survminer	Drawing Survival Curves using 'ggplot2'
survMisc	Miscellaneous Functions for Survival Data
survPresmooth	Presmoothed Estimation in Survival Analysis
SurvRank	Rank Based Survival Modelling
survrec	Survival analysis for recurrent event data
SurvRegCensCov	Weibull Regression for a Right-Censored Endpoint with Interval-Censored Covariate
survRM2	Comparing Restricted Mean Survival Time
survsim	Simulation of Simple and Complex Survival Data
survSNP	Power Calculations for SNP Studies with Censored Outcomes
sValues	Measures of the Sturdiness of Regression Coefficients
svapls	Surrogate variable analysis using partial least squares in a gene expression study
svcm	2d and 3d Space-Varying Coefficient Models
svd	Interfaces to Various State-of-Art SVD and Eigensolvers
svDialogs	SciViews GUI API - Dialog boxes
svDialogstcltk	SciViews GUI API - Dialog boxes using Tcl/Tk
svdvis	Singular Value Decomposition Visualization
svdvisual	SVD visualization tools
svglite	An 'SVG' Graphics Device
svgPanZoom	R 'Htmlwidget' to Add Pan and Zoom to Almost any R Graphic
svGUI	SciViews GUI API - Functions to manage GUIs
svgViewR	3D Animated Interactive Visualizations Using SVG
svHttp	SciViews GUI API - R HTTP server
svIDE	SciViews GUI API - IDE and code editor functions
svKomodo	SciViews GUI API - Functions to interface with Komodo Edit/IDE
svmadmm	Linear/Nonlinear SVM Classification Solver Based on ADMM and IADMM Algorithms
svMisc	SciViews GUI API - Miscellaneous functions
SVMMaj	SVMMaj algorithm
SVMMatch	Causal Effect Estimation and Diagnostics with Support Vector Machines
svmpath	svmpath: the SVM Path algorithm
svs	Tools for Semantic Vector Spaces
svSocket	SciViews GUI API - R Socket Server
svSweave	SciViews GUI API - Sweave functions
svTools	SciViews GUI API - Tools (wrapper for packages tools and codetools)
svUnit	SciViews GUI API - Unit testing
svWidgets	SciViews GUI API - Widgets & Windows
SvyNom	Nomograms for Right-Censored Outcomes from Survey Designs
svyPVpack	A package for complex surveys including plausible values
swamp	Visualization, analysis and adjustment of high-dimensional data in respect to sample annotations
SwarmSVM	Ensemble Learning Algorithms Based on Support Vector Machines
SWATmodel	A multi-OS implementation of the TAMU SWAT model
SweaveListingUtils	Utilities for Sweave together with TeX listings package
sweidnumbr	Handling of Swedish Identity Numbers
swfscMisc	Miscellaneous Functions for Southwest Fisheries Science Center
swirl	Learn R, in R
swirlify	A Toolbox for Writing 'swirl' Courses
SwissAir	Air Quality Data of Switzerland for one year in 30 min Resolution
switchnpreg	Switching nonparametric regression models for a single curve and functional data
switchr	Installing, Managing, and Switching Between Distinct Sets of Installed Packages
switchrGist	Publish Package Manifests to GitHub Gists
SWMPr	Retrieving, Organizing, and Analyzing Estuary Monitoring Data
sybil	Efficient Constrained Based Modelling in R
sybilccFBA	Cost Constrained FLux Balance Analysis: MetabOlic Modeling with ENzyme kineTics (MOMENT)
sybilcycleFreeFlux	cycle-Free Flux balance analysis: Efficient removal of thermodynamically infeasible cycles from metabolic flux distributions
sybilDynFBA	Dynamic FBA : dynamic flux balance analysis
sybilEFBA	Using Gene Expression Data to Improve Flux Balance Analysis Predictions
sybilSBML	SBML Integration in Package Sybil
symbolicDA	Analysis of Symbolic Data
symbols	Symbol plots
symmoments	Symbolic central and noncentral moments of the multivariate normal distribution
synbreed	Framework for the Analysis of Genomic Prediction Data using R
synbreedData	Data for the Synbreed Package
synchronicity	Boost Mutex Functionality in R
synchrony	Methods for computing spatial, temporal, and spatiotemporal statistics
SynchWave	Synchrosqueezed Wavelet Transform
SyncMove	Subsample Temporal Data to Synchronal Events and Compute the MCI
SYNCSA	SYNCSA - Analysis of functional and phylogenetic patterns in metacommunities
SynergizeR	Interface to The Synergizer service for translating between sets of biological identifiers
SyNet	Inference and Analysis of Sympatry Networks
synlik	Synthetic Likelihood methods for intractable likelihoods
synRNASeqNet	Synthetic RNA-Seq Network Generation and Mutual Information Estimates
Synth	Synthetic Control Group Method for Comparative Case Studies
synthpop	Generating Synthetic Versions of Sensitive Microdata for Statistical Disclosure Control
sysfonts	Loading System Fonts into R
systemfit	Estimating Systems of Simultaneous Equations
systemicrisk	A Toolbox for Systemic Risk
syuzhet	Extracts Sentiment and Sentiment-Derived Plot Arcs from Text
0	
tab	Functions for Creating Summary Tables for Statistical Reports
taber	Split and Recombine Your Data
tablaxlsx	Write Formatted Tables in Excel Workbooks
Table1Heatmap	Table 1 Heatmap
table1xls	Produces summary tables and exports them to multi-tab spreadsheet format (.xls or .xlsx)
TableMonster	Table Monster
tableone	Create "Table 1" to Describe Baseline Characteristics
tableplot	Represents tables as semi-graphic displays
tables	Formula-driven table generation
TableToLongForm	TableToLongForm
tabplot	Tableplot, a Visualization of Large Datasets
tabplotd3	Tabplotd3, interactive inspection of large data
tabuSearch	R based tabu search algorithm
tadaatoolbox	Helpers for Data Analysis and Presentation Focused on Undergrad Psychology
tagcloud	Tag Clouds
tailDepFun	Minimum Distance Estimation of Tail Dependence Models
tailloss	Estimate the Probability in the Upper Tail of the Aggregate Loss Distribution
TAM	Test Analysis Modules
TANOVA	Time Course Analysis of Variance for Microarray
TaoTeProgramming	Illustrations from Tao Te Programming
TapeR	Flexible Tree Taper Curves Based on Semiparametric Mixed Models
TAQMNGR	Manage Tick-by-Tick Transaction Data
Tariff	Replicate Tariff Method for Verbal Autopsy
taRifx	Collection of utility and convenience functions
taRifx.geo	Collection of various spatial functions
tau	Text Analysis Utilities
TauP.R	Earthquake Traveltime Calculations for 1-D Earth Models
TauStar	Efficient Computation and Testing of the Bergsma-Dassios Sign Covariance
tawny	Provides various portfolio optimization strategies including random matrix theory and shrinkage estimators
tawny.types	Common types for tawny
taxize	Taxonomic Information from Around the Web
Taxonstand	Taxonomic Standardization of Plant Species Names
tbart	Teitz and Bart's p-Median Algorithm
tbdiag	Functions for tuberculosis diagnostics research
TBEST	Tree Branches Evaluated Statistically for Tightness
TBSSurvival	Survival Analysis using a Transform-Both-Sides model
TCGA2STAT	Simple TCGA Data Access for Integrated Statistical Analysis in R
TcGSA	Time-Course Gene Set Analysis
tcltk2	Tcl/Tk Additions
tclust	Robust Trimmed Clustering
tcR	Advanced Data Analysis of Immune Receptor Repertoires
TDA	Statistical Tools for Topological Data Analysis
TDAmapper	Analyze High-Dimensional Data Using Discrete Morse Theory
TDboost	A Boosted Tweedie Compound Poisson Model
TDCor	Gene Network Inference from Time-Series Transcriptomic Data
TDD	Time-Domain Deconvolution of Seismometer Response
TDMR	Tuned Data Mining in R
tdr	Target Diagram
tdthap	TDT tests for extended haplotypes
TeachingDemos	Demonstrations for Teaching and Learning
TeachingSampling	Selection of Samples and Parameter Estimation in Finite Population
TeachNet	Fits neural networks to learn about back propagation
TED	Turbulence Time Series Event Detection and Classification
TEEReg	Trimmed Elemental Estimation for Linear Models
teigen	Model-Based Clustering and Classification with the Multivariate t Distribution
telegram	R Wrapper Around the Telegram Bot API
tempcyclesdata	Climate Data from Wang and Dillon
tempdisagg	Methods for Temporal Disaggregation and Interpolation of Time Series
tensor	Tensor product of arrays
tensorA	Advanced tensors arithmetic with named indices
tensr	Covariance Inference and Decompositions for Tensor Datasets
TEQR	Target Equivalence Range Design
TERAplusB	Test for A+B Traditional Escalation Rule
tergm	Fit, Simulate and Diagnose Models for Network Evolution Based on Exponential-Family Random Graph Models
termstrc	Zero-coupon Yield Curve Estimation
ternvis	Visualisation, verification and calibration of ternary probabilistic forecasts
TESS	Diversification Rate Estimation and Fast Simulation of Reconstructed Phylogenetic Trees under Tree-Wide Time-Heterogeneous Birth-Death Processes Including Mass-Extinction Events
tester	Tests and checks characteristics of R objects
TestingSimilarity	Bootstrap Test for Similarity of Dose Response Curves Concerning the Maximum Absolute Deviation
testit	A Simple Package for Testing R Packages
TestScorer	GUI for Entering Test Items and Obtaining Raw and Transformed Scores
TestSurvRec	Statistical tests to compare two survival curves with recurrent events
testthat	Unit Testing for R
testthatsomemore	A mocking, stubbing, and file testing framework extending testthat
texmex	Statistical modelling of extreme values
texmexseq	Treatment Effect eXplorer for Microbial Ecology eXperiments (using Sequence Counts)
TExPosition	Two-table ExPosition
texreg	Conversion of R Regression Output to LaTeX or HTML Tables
text2vec	Fast and Modern Text Mining Framework - Vectorization and Word Embeddings
textcat	N-Gram Based Text Categorization
textir	Inverse Regression for Text Analysis
textmineR	Functions for Text Mining and Topic Modeling
textometry	Textual Data Analysis Package used by the TXM Software
TextoMineR	Textual Statistics
textreg	n-Gram Text Regression, aka Concise Comparative Summarization
textreuse	Detect Text Reuse and Document Similarity
TFDEA	Technology Forecasting using DEA (Data Envelopment Analysis)
tfer	Forensic Glass Transfer Probabilities
TFMPvalue	Efficient and Accurate P-Value Computation for Position Weight Matrices
tfplot	Time Frame User Utilities
tframe	Time Frame Coding Kernel
tframePlus	Time Frame Coding Kernel Extensions
TFX	R API to TrueFX(tm)
tgcd	Thermoluminescence Glow Curve Deconvolution
tggd	The Standard Distribution Functions for the Truncated Generalised Gamma Distribution
tglm	Binary Regressions under Independent Student-t Priors
tgp	Bayesian Treed Gaussian Process Models
tgram	Functions to compute and plot tracheidograms
TH.data	TH's Data Archive
Thermimage	Functions for Handling Thermal Images
thermocouple	Temperature Measurement with Thermocouples, RTD and IC Sensors
thgenetics	Genetic Rare Variants Tests
Thinknum	Thinknum Data Connection
ThreeArmedTrials	Design and Analysis of Clinical Non-inferiority or Superiority Trials with Active and Placebo Control
threeboost	Thresholded variable selection and prediction based on estimating equations
ThreeGroups	ML Estimator for Baseline-Placebo-Treatment (Three-Group) Experiments
threejs	Interactive 3D Scatter Plots and Globes
ThreeWay	Three-Way Component Analysis
threewords	Represent Precise Coordinates in Three Words
threg	Threshold Regression
thregI	Threshold Regression for Interval-Censored Data with Cure-Rate or without Cure-Rate Model
ThresholdROC	Threshold Estimation
thsls	Three-Stage Least Squares Estimation for Systems of Simultaneous Equations
tibble	Simple Data Frames
tibbrConnector	R interface to tibbr
TickExec	Execution Functions for Tick Data Back Test
tictoc	Functions for timing R scripts, as well as implementations of Stack and List structures
TiddlyWikiR	Create dynamic reports using a TiddlyWiki template
TideCurves	Analysis and Prediction of Tides
TideHarmonics	Harmonic Analysis of Tides
Tides	Quasi-Periodic Time Series Characteristics
TideTables	Tide Analysis and Prediction of Predominantly Semi-Diurnal Tides
tidyjson	A Grammar for Turning 'JSON' into Tidy Tables
tidyr	Easily Tidy Data with 'spread()' and 'gather()' Functions
tiff	Read and write TIFF images
tiger	TIme series of Grouped ERrors
tigerstats	R Functions for Elementary Statistics
tightClust	Tight Clustering
tigris	Load Census TIGER/Line Shapefiles into R
tikzDevice	R Graphics Output in LaTeX Format
tileHMM	Hidden Markov Models for ChIP-on-Chip Analysis
TilePlot	Characterization of functional genes in complex microbial communities using tiling DNA microarrays
timeDate	Rmetrics - Chronological and Calendar Objects
timedelay	Time Delay Estimation for Stochastic Time Series of Gravitationally Lensed Quasars
timeit	Easy profiling of R functions
timeline	Timelines for a Grammar of Graphics
TimeMachine	Time Machine
timeordered	Time-ordered and time-aggregated network analyses
TimeProjection	Time Projections
timereg	Flexible Regression Models for Survival Data
timeROC	Time-Dependent ROC Curve and AUC for Censored Survival Data
timesboot	Bootstrap computations for time series objects
timeSeq	Detecting Differentially Expressed Genes in Time Course RNA-Seq Data
timeSeries	Rmetrics - Financial Time Series Objects
timeseriesdb	Manage Time Series with R and PostgreSQL
timetools	Seasonal/Sequential (Instants/Durations, Even or not) Time Series
timetree	Interface to the TimeTree of Life Webpage
TimeWarp	Date Calculations and Manipulation
timma	Target Inhibition Interaction using Maximization and Minimization Averaging
TIMP	Fitting Separable Nonlinear Models in Spectroscopy and Microscopy
timsac	TIMe Series Analysis and Control package
Tinflex	Tinflex - Universal Non-Uniform Random Number Generator
TInPosition	Inference tests for TExPosition
TipDatingBeast	Using Tip Dates with Phylogenetic Trees in BEAST
tipom	Automated measure-based classification for flint tools
tis	Time Indexes and Time Indexed Series
titan	Titration analysis for mass spectrometry data
TITAN2	Threshold Indicator Taxa Analysis
titanic	Titanic Passenger Survival Data Set
titrationCurves	Acid/Base, Complexation, Redox, and Precipitation Titration Curves
TKF	Pairwise Distance Estimation with TKF91 and TKF92 Model
tkrgl	TK widget tools for rgl package
tkrplot	TK Rplot
TLBC	Two-Level Behavior Classification
TLdating	Tools for Thermoluminescences Dating
tlemix	Trimmed Maximum Likelihood Estimation
tlm	Effects under Linear, Logistic and Poisson Regression Models with Transformed Variables
tlmec	Linear Student-t Mixed-Effects Models with Censored Data
tlnise	Two-level normal independent sampling estimation
tm	Text Mining Package
tmap	Thematic Maps
TMB	Template Model Builder: A General Random Effect Tool Inspired by ADMB
TMDb	Access to TMDb API - Apiary
tmg	Truncated Multivariate Gaussian Sampling
Tmisc	Turner Miscellaneous
tmle	Targeted Maximum Likelihood Estimation
tmlenet	Targeted Maximum Likelihood Estimation for Network Data
tmle.npvi	Targeted Learning of a NP Importance of a Continuous Exposure
tmod	Transcriptional Module Analysis
tm.plugin.alceste	Import texts from files in the Alceste format using the tm text mining framework
tm.plugin.dc	Text Mining Distributed Corpus Plug-In
tm.plugin.europresse	Import Articles from 'Europresse' Using the 'tm' Text Mining Framework
tm.plugin.factiva	Import articles from Factiva using the tm text mining framework
tm.plugin.lexisnexis	Import Articles from LexisNexis Using the tm Text Mining Framework
tm.plugin.mail	Text Mining E-Mail Plug-In
tm.plugin.webmining	Retrieve Structured, Textual Data from Various Web Sources
tmpm	Trauma Mortality Prediction Model
tmvnsim	Truncated Multivariate Normal Simulation
tmvtnorm	Truncated Multivariate Normal and Student t Distribution
tnam	Temporal Network Autocorrelation Models (TNAM)
tnet	Software for Analysis of Weighted, Two-Mode, and Longitudinal Networks
toaster	Big Data in-Database Analytics with Teradata Aster
TOC	Total Operating Characteristic Curve and ROC Curve
tolBasis	Fundamental Definitions and Utilities of the Time Oriented Language (TOL)
tolerance	Functions for Calculating Tolerance Intervals
toOrdinal	Function for Converting Cardinal to Ordinal Numbers by Adding a Language Specific Ordinal Indicator to the Number
topicmodels	Topic Models
TopKLists	Inference, Aggregation and Visualization for Top-K Ranked Lists
topmodel	Implementation of the hydrological model TOPMODEL in R
topologyGSA	Gene Set Analysis Exploiting Pathway Topology
topsis	TOPSIS method for multiple-criteria decision making (MCDM)
tosls	Instrumental Variables Two Stage Least Squares estimation
TotalCopheneticIndex	Total Cophenetic Index
tourr	Implement Tour Methods in R Code
tourrGui	A Tour GUI using gWidgets
toxtestD	Experimental design for binary toxicity tests
tpe	Tree preserving embedding
TP.idm	Estimation of Transition Probabilities for the Illness-Death Model
TPmsm	Estimation of Transition Probabilities in Multistate Models
tpr	Temporal Process Regression
TR8	A Tool for Downloading Functional Traits Data for Plant Species
tracheideR	Standardize Tracheidograms
track	Store Objects on Disk Automatically
trackeR	Infrastructure for Running and Cycling Data from GPS-Enabled Tracking Devices
TrackReconstruction	Reconstruct animal tracks from magnetometer, accelerometer, depth and optional speed data
tractor.base	A package for reading, manipulating and visualising magnetic resonance images
TRADER	Tree Ring Analysis of Disturbance Events in R
traitr	An interface for creating GUIs modeled in part after traits UI module for python
traits	Species Trait Data from Around the Web
Traitspace	A Predictive Model for Trait Based Community Assembly of Plant Species
traj	Trajectory Analysis
trajectories	Classes and Methods for Trajectory Data
TraMineR	Trajectory Miner: a Toolbox for Exploring and Rendering Sequences
TraMineRextras	TraMineR Extension
TRAMPR	TRFLP Analysis and Matching Package for R
transcribeR	Automated Transcription of Audio Files Through the HP IDOL API
TransferEntropy	The Transfer Entropy Package
translate	Bindings for the Google Translate API v2
translateR	Bindings for the Google and Microsoft Translation APIs
translateSPSS2R	Toolset for Translating SPSS-Syntax to R-Code
translation.ko	R Manuals Literally Translated in Korean
TransModel	Fit Linear Transformation Models for Censored Data
TransP	Implementation of Transportation Problem Algorithms
transport	Optimal Transport in Various Forms
trapezoid	The Trapezoidal Distribution
TRD	Transmission Ratio Distortion
TreatmentSelection	Evaluate Treatment Selection Biomarkers
treatSens	Sensitivity Analysis for Causal Inference
tree	Classification and Regression Trees
treebase	Discovery, Access and Manipulation of 'TreeBASE' Phylogenies
treeclim	Numerical Calibration of Proxy-Climate Relationships
treeClust	Cluster Distances Through Trees
treecm	Centre of Mass Assessment and Consolidation of Trees
treelet	An Adaptive Multi-Scale Basis for High-Dimensional, Sparse and Unordered Data
treemap	Treemap Visualization
TreePar	Estimating birth and death rates based on phylogenies
treeperm	Exact and Asymptotic K Sample Permutation Test
treeplyr	'dplyr' Functionality for Matched Tree and Data Objects
treescape	Statistical Exploration of Landscapes of Phylogenetic Trees
TreeSim	Simulating Phylogenetic Trees
TreeSimGM	Simulating Phylogenetic Trees under a General Model
treethresh	Methods for Tree-based Local Adaptive Thresholding
trelliscope	Create and Navigate Large Multi-Panel Visual Displays
trend	Non-Parametric Trend Tests and Change-Point Detection
TrialSize	R functions in Chapter 3,4,6,7,9,10,11,12,14,15
triangle	Provides the Standard Distribution Functions for the Triangle Distribution
trib	Analysing and Visualizing Tribology Measurements
trifield	Some basic facilities for ternary fields and plots
TriMatch	Propensity Score Matching of Non-Binary Treatments
trimcluster	Cluster analysis with trimming
trimr	An Implementation of Common Response Time Trimming Methods
trimTrees	Trimmed opinion pools of trees in a random forest
trioGxE	A data smoothing approach to explore and test gene-environment interaction in case-parent trio data
trip	Tools for the Analysis of Animal Track Data
tripack	Triangulation of Irregularly Spaced Data
tripEstimation	Metropolis Sampler and Supporting Functions for Estimating Animal Movement from Archival Tags and Satellite Fixes
TripleR	Social Relation Model (SRM) Analyses for Single or Multiple Groups
TROM	Transcriptome Overlap Measure
trotter	Pseudo-Vectors Containing All Permutations, Combinations and Subsets of Objects Taken from a Vector
TRSbook	Functions and Datasets to Accompany the Book "The R Software: Fundamentals of Programming and Statistical Analysis"
trueskill	Implementation the TrueSkill algorithm in R
TruncatedNormal	Truncated Multivariate Normal
truncdist	Truncated Random Variables
truncgof	GoF tests allowing for left truncated data
truncnorm	Truncated normal distribution
truncreg	Truncated Gaussian Regression Models
truncSP	Semi-parametric estimators of truncated regression models
trust	Trust Region Optimization
trustOptim	Trust Region Optimization for Nonlinear Functions with Sparse Hessians
TSA	Time Series Analysis
tsallisqexp	Tsallis q-Exp Distribution
tsbridge	Calculate normalising constants for Bayesian time series models
tsBSS	Tools for Blind Source Separation for Time Series
tsbugs	Create time series BUGS models
tsc	Likelihood-ratio Tests for Two-Sample Comparisons
TSclust	Time Series Clustering Utilities
TScompare	'TSdbi' Database Comparison
tscount	Analysis of Count Time Series
TSdata	'TSdbi' Illustration
TSdbi	Time Series Database Interface
TSdist	Distance Measures for Time Series Data
tsDyn	Nonlinear Time Series Models with Regime Switching
tseries	Time Series Analysis and Computational Finance
tseriesChaos	Analysis of nonlinear time series
tseriesEntropy	Entropy Based Analysis and Tests for Time Series
tsfa	Time Series Factor Analysis
TSfame	'TSdbi' Extensions for Fame
TSHRC	Two Stage Hazard Rate Comparison
tsintermittent	Intermittent Time Series Forecasting
TSMining	Mining Univariate and Multivariate Motifs in Time-Series Data
TSmisc	'TSdbi' Extensions to Wrap Miscellaneous Data Sources
tsModel	Time Series Modeling for Air Pollution and Health
TSMySQL	'TSdbi' Extensions for 'MySQL'
tsna	Tools for Temporal Social Network Analysis
tsne	T-distributed Stochastic Neighbor Embedding for R (t-SNE)
TSodbc	'TSdbi' Extensions for ODBC
tsoutliers	Detection of Outliers in Time Series
TSP	Traveling Salesperson Problem (TSP)
Tsphere	Transposable Sphering for Large-Scale Inference with Correlated Data
tsPI	Improved Prediction Intervals for ARIMA Processes and Structural Time Series
tspmeta	Instance Feature Calculation and Evolutionary Instance Generation for the Traveling Salesman Problem
TSPostgreSQL	'TSdbi' Extensions for 'PostgreSQL'
TSPred	Functions for Baseline-Based Time Series Prediction
TSsdmx	'TSdbi' Extension to Connect with 'SDMX'
TSsql	Generic SQL Helper Functions for 'TSdbi' SQL Plugins
TSSQLite	'TSdbi' Extensions for 'SQLite'
TSTr	Ternary Search Tree for Auto-Completion and Spell Checking
TSTutorial	Fitting and Predict Time Series Interactive Laboratory
TTAinterfaceTrendAnalysis	Temporal Trend Analysis Graphical Interface
tth	TeX to HTML/MathML Translators tth/ttm
TTmoment	Sampling and Calculating the First and Second Moments for the Doubly Truncated Multivariate t Distribution
TTR	Technical Trading Rules
TTS	Master Curve Estimates Corresponding to Time-Temperature Superposition
ttScreening	Genome-wide DNA methylation sites screening by use of training and testing samples
tttplot	Time to Target Plot
ttutils	Utility functions
ttwa	Travel To Work Area
tuber	Client for the YouTube API
tufte	Tufte's Styles for R Markdown Documents
tufterhandout	Tufte-style html document format for rmarkdown
TukeyC	Conventional Tukey Test
tumblR	Access to Tumblr v2 API
tumgr	Tumor Growth Rate Analysis
TunePareto	Multi-objective parameter tuning for classifiers
tuneR	Analysis of music and speech
tuple	Find every match, or orphan, duplicate, triplicate, or other replicated values
turboEM	A Suite of Convergence Acceleration Schemes for EM, MM and other fixed-point algorithms
turfR	TURF Analysis for R
turner	Turn vectors and lists of vectors into indexed structures
TurtleGraphics	Turtle graphics in R
tutorial	Convert R Markdown Files to DataCamp Light HTML Files
TUWmodel	Lumped Hydrological Model for Education Purposes
tvd	Total Variation Denoising
tvm	Time Value of Money Functions
twang	Toolkit for Weighting and Analysis of Nonequivalent Groups
tweedie	Tweedie exponential family models
tweenr	Interpolate Data for Smooth Animations
tweet2r	Twitter Collector and Export to 'SQLite', 'postGIS' and 'GIS' Format
twiddler	Interactive manipulation of R expressions
twitteR	R Based Twitter Client
TwoCop	Nonparametric test of equality between two copulas
TwoPhaseInd	Estimate Gene-Treatment Interaction Exploiting Randomization
twoStageGwasPower	Compute thresholds and power for two-stage gwas
twostageTE	Two-Stage Threshold Estimation
TwoStepCLogit	Conditional Logistic Regression: A Two-Step Estimation Method
txtplot	Text based plots
0	
UBCRM	Functions to Simulate and Conduct Dose-Escalation Phase I Studies
ucbthesis	UC Berkeley graduate division thesis template
ucminf	General-purpose unconstrained non-linear optimization
udunits2	Udunits-2 Bindings for R
ukgasapi	API for UK Gas Market Information
Ultimixt	Bayesian Analysis of a Non-Informative Parametrization for Gaussian Mixture Distributions
ump	Uniformly Most Powerful Tests
umx	Structural Equation Modelling in R with OpenMx
unbalanced	Racing for Unbalanced Methods Selection
unbalhaar	Function estimation via Unbalanced Haar wavelets
UncerIn2	Implements Models of Uncertainty into the Interpolation Functions
UNF	Tools for Creating Universal Numeric Fingerprints for Data
unfoldr	Stereological Unfolding for Spheroidal Particles
uniah	Unimodal Additive Hazards Model
Unicode	Unicode Data and Utilities
uniCox	Univarate shrinkage prediction in the Cox model
uniftest	Tests for Uniformity
uniqtag	Abbreviate Strings to Short, Unique Identifiers
uniqueAtomMat	Finding Unique or Duplicated Rows or Columns for Atomic Matrices
uniReg	Unimodal penalized spline regression using B-splines
unitedR	Assessment and Evaluation of Formations in United
unittest	TAP-Compliant Unit Testing
unmarked	Models for Data from Unmarked Animals
untb	ecological drift under the UNTB
upclass	Updated Classification Methods using Unlabeled Data
uplift	Uplift Modeling
UPMASK	Unsupervised Photometric Membership Assignment in Stellar Clusters
UpSetR	A More Scalable Alternative to Venn and Euler Diagrams for Visualizing Intersecting Sets
uptimeRobot	Access the UptimeRobot Ping API
urca	Unit Root and Cointegration Tests for Time Series Data
urlshorteneR	R Wrapper for the 'Bit.ly', 'Goo.gl' and 'Is.gd' URL Shortening Services
urltools	Vectorised Tools for URL Handling and Parsing
uroot	Unit Root Tests for Seasonal Time Series
USAboundaries	Historical and Contemporary Boundaries of the United States of America
UScancer	Create US cancer datasets from SEER, IARC, and US Census data
UScensus2000cdp	US Census 2000 Designated Places Shapefiles and Additional Demographic Data
UScensus2000tract	US Census 2000 Tract Level Shapefiles and Additional Demographic Data
UScensus2010	US Census 2010 Suite of R Packages
usdm	Uncertainty Analysis for Species Distribution Models
useful	A Collection of Handy, Useful Functions
userfriendlyscience	Quantitative Analysis Made Accessible
UsingR	Data Sets, Etc. for the Text "Using R for Introductory Statistics", Second Edition
uskewFactors	Model-based clustering via mixtures of unrestricted skew-t factor analyzer models
usl	Analyze System Scalability with the Universal Scalability Law
ustyc	Fetch US Treasury yield curve data
utility	Construct, Evaluate and Plot Value and Utility Functions
uuid	Tools for generating and handling of UUIDs
UWHAM	Unbinned weighted histogram analysis method (UWHAM)
0	
V8	Embedded JavaScript Engine
vacem	Vaccination Activities Coverage Estimation Model
validate	Data Validation Infrastructure
validateRS	One-Sided Multivariate Testing Procedures for Rating Systems
valottery	Results from the Virginia Lottery Draw Games
varbvs	Variational inference for Bayesian variable selection
varComp	Variance Component Models
vardiag	Variogram Diagnostics
vardpoor	Variance Estimation for Sample Surveys by the Ultimate Cluster Method
VaRES	Computes value at risk and expected shortfall for over 100 parametric distributions
VAR.etp	VAR modelling: estimation, testing, and prediction
VarfromPDB	Capture the Genes and Variants Related to a Genetic Disease from Public Databases
varhandle	Functions for Robust Variable Handling
VariABEL	Testing of genotypic variance heterogeneity to detect potentially interacting SNP
variables	Variable Descriptions
varian	Variability Analysis in R
VarianceGamma	The Variance Gamma Distribution
vars	VAR Modelling
VARSEDIG	An Algorithm for Morphometric Characters Selection and Statistical Validation in Morphological Taxonomy
VarSelLCM	Variable Selection for Model-Based Clustering using the Integrated Complete-Data Likelihood of a Latent Class Model
varSelRF	Variable Selection using Random Forests
VARsignR	Sign Restrictions, Bayesian, Vector Autoregression Models
VarSwapPrice	Pricing a variance swap on an equity index
vartors	Transform Definition of Variables to R Scripts
vbdm	Variational Bayes Discrete Mixture Model
VBLPCM	Variational Bayes Latent Position Cluster Model for Networks
VBmix	Variational Bayesian Mixture Models
vbsr	Variational Bayes Spike Regression Regularized Linear Models
VCA	Variance Component Analysis
vcd	Visualizing Categorical Data
vcdExtra	'vcd' Extensions and Additions
vcfR	Tools for Working with Variant Call Format ('VCF') Files
vcrpart	Tree-Based Varying Coefficient Regression for Generalized Linear and Ordinal Mixed Models
VDA	VDA
vdg	Variance Dispersion Graphs and Fraction of Design Space Plots
Vdgraph	Variance dispersion graphs and Fraction of design space plots for response surface designs
VdgRsm	Plots of Scaled Prediction Variances for Response Surface Designs
vdmR	Visual Data Mining Tools for R
vec2dtransf	2D Cartesian Coordinate Transformation
vecsets	like base::sets tools but keeps duplicate elements
VecStatGraphs2D	Vector analysis using graphical and analytical methods in 2D
VecStatGraphs3D	Vector analysis using graphical and analytical methods in 3D
vegalite	Tools to Encode Visualizations with the 'Grammar of Graphics'-Like 'Vega-Lite' 'Spec'
vegan	Community Ecology Package
vegan3d	Static and Dynamic 3D Plots for the 'vegan' Package
vegclust	Fuzzy clustering of vegetation data
vegdata	Access Vegetation Databases and Treat Taxonomy
vegetarian	Jost Diversity Measures for Community Data
venn	Draw Venn Diagrams
VennDiagram	Generate High-Resolution Venn and Euler Plots
venneuler	Venn and Euler Diagrams
verification	Weather Forecast Verification Utilities
versions	Query and Install Specific Versions of Packages on CRAN
vertexenum	Vertex Enumeration of Polytopes
VertexSimilarity	Creates Vertex Similarity Matrix for an Undirected Graph
vetools	Tools for Venezuelan Environmental Data
VGAM	Vector Generalized Linear and Additive Models
VGAMdata	Data Supporting the 'VGAM' Package
VHDClassification	Discrimination/Classification in very high dimension with linear and quadratic rules
VideoComparison	Video Comparison Tool
VIF	VIF Regression: A Fast Regression Algorithm For Large Data
VIFCP	Detecting Change-Points via VIFCP Method
VIGoR	Variational Bayesian Inference for Genome-Wide Regression
VIM	Visualization and Imputation of Missing Values
VIMGUI	Visualization and Imputation of Missing Values
VineCopula	Statistical Inference of Vine Copulas
vines	Multivariate Dependence Modeling with Vines
violinmplot	Combination of violin plot with mean and standard deviation
vioplot	Violin plot
viopoints	1-D Scatter Plots with Jitter Using Kernel Density Estimates
vipor	Plot Categorical Data Using Quasirandom Noise and Density Estimates
viridis	Default Color Maps from 'matplotlib'
viridisLite	Default Color Maps from 'matplotlib' (Lite Version)
virtualspecies	Generation of Virtual Species Distributions
ViSiElse	A Visual Tool for Behaviour Analysis
visNetwork	Network Visualization using 'vis.js' Library
visreg	Visualization of Regression Models
visualFields	Statistical Methods for Visual Fields
visualize	Graph Probability Distributions with User Supplied Parameters and Stats
VisuClust	Visualisation of Clusters in Multivariate Data
vita	Variable Importance Testing Approaches
vitality	Fitting Routines for the Vitality Family of Mortality Models
VizOR	Graphical visualization tools for complex observational data with focus on health sciences (VizOR)
VLF	Frequency Matrix Approach for Assessing Very Low Frequency Variants in Sequence Records
VLMC	Variable Length Markov Chains ('VLMC') Models
vmsbase	GUI Tools to Process, Analyze and Plot Fisheries Data
VNM	Using V-algorithm and Newton-Raphson Method to Obtain Multiple-objective Optimal Design
Voss	Generic Voss algorithm (random sequential additions)
vottrans	Voter Transition Analysis
vowels	Vowel Manipulation, Normalization, and Plotting
vows	Voxelwise Semiparametrics
VoxR	Metrics extraction of trees from T-LiDAR data
VPdtw	Variable Penalty Dynamic Time Warping
vqtl	Genome Scans to Accommodate and Target Genetic and Non-Genetic Effects on Trait Variance
vrcp	Change Point Estimation for Regression with Varying Segments and Heteroscedastic Variances
vrmlgen	Generate 3D visualizations for data exploration on the web
vrtest	Variance Ratio tests and other tests for Martingale Difference Hypothesis
vscc	Variable selection for clustering and classification
VSE	Variant Set Enrichment
VSURF	Variable Selection Using Random Forests
VTrack	A Collection of Tools for the Analysis of Remote Acoustic Telemetry Data
vtreat	Simple Variable Treatment
vudc	Visualization of Univariate Data for Comparison
vwr	Useful functions for visual word recognition research
0	
W2CWM2C	A Graphical Tool for Wavelet (Cross) Correlation and Wavelet Multiple (Cross) Correlation Analysis
W3CMarkupValidator	R Interface to W3C Markup Validation Services
WACS	Multivariate Weather-State Approach Conditionally Skew-Normal Generator
waffect	A package to simulate constrained phenotypes under a disease model H1
waffle	Create Waffle Chart Visualizations in R
wahc	Autocorrelation and Heteroskedasticity Correction in Fixed Effect Panel Data Model
wakefield	Generate Random Data Sets
walkr	Random Walks in the Intersection of Hyperplanes and the N-Simplex
walkscoreAPI	Walk Score and Transit Score API
warbleR	Streamline Bioacoustic Analysis
WARN	Weaning Age Reconstruction with Nitrogen Isotope Analysis
wasim	Visualisation and analysis of output files of the hydrological model WASIM
water	Actual Evapotranspiration with Energy Balance Models
waterData	An R Package for Retrieval, Analysis, and Anomaly Calculation of Daily Hydrologic Time Series Data
WaterML	Fetch and Analyze Data from 'WaterML' and 'WaterOneFlow' Web Services
Watersheds	Spatial Watershed Aggregation and Spatial Drainage Network Analysis
Wats	Wrap Around Time Series Graphics
waveband	Computes credible intervals for Bayesian wavelet shrinkage
waved	Wavelet Deconvolution
WaveletComp	Computational Wavelet Analysis
wavelets	A package of functions for computing wavelet filters, wavelet transforms and multiresolution analyses
wavemulcor	Wavelet routine for multiple correlation
WaverR	Data Estimation using Weighted Averages of Multiple Regressions
waveslim	Basic wavelet routines for one-, two- and three-dimensional signal processing
wavethresh	Wavelets statistics and transforms
wBoot	Bootstrap Methods
wbs	Wild Binary Segmentation for Multiple Change-Point Detection
wbsts	Multiple Change-Point Detection for Nonstationary Time Series
wccsom	SOM Networks for Comparing Patterns with Peak Shifts
WCE	Weighted Cumulative Exposure Models
WCQ	Detection of QTL effects in a small mapping population
WDI	World Development Indicators (World Bank)
weatherData	Get Weather Data from the Web
weathermetrics	Functions to Convert Between Weather Metrics
weatherr	Tools for Handling and Scrapping Instant Weather Forecast Feeds
webchem	Chemical Information from the Web
webp	A New Format for Lossless and Lossy Image Compression
webreadr	Tools for Reading Formatted Access Log Files
webshot	Take Screenshots of Web Pages
webuse	Import Stata 'webuse' Datasets
webutils	Utility Functions for Web Applications
webvis	Create graphics for the web from R
wec	Weighted Effect Coding
WeightedCluster	Clustering of Weighted Data
Weighted.Desc.Stat	Weighted Descriptive Statistics
WeightedPortTest	Weighted Portmanteau Tests for Time Series Goodness-of-fit
weightedScores	Weighted Scores Method for Regression Models with Dependent Data
weights	Weighting and Weighted Statistics
weightTAPSPACK	Weight TAPS Data
weirs	A Hydraulics Package to Compute Open-Channel Flow over Weirs
wellknown	Convert Between 'WKT' and 'GeoJSON'
wesanderson	A Wes Anderson Palette Generator
wfe	Weighted Linear Fixed Effects Regression Models for Causal Inference
wfg	Weighted Fast Greedy Algorithm
wgaim	Whole Genome Average Interval Mapping for QTL Detection using Mixed Models
WGCNA	Weighted Correlation Network Analysis
wgsea	Wilcoxon based gene set enrichment analysis
WhatIf	WhatIf: Software for Evaluating Counterfactuals
whisker	{{mustache}} for R, logicless templating
WhiteStripe	White Matter Normalization for Magnetic Resonance Images using Whitestripe
WHO	R Client for the World Health Organization API
whoami	Username, Full Name, Email Address, 'GitHub' Username of the Current User
whoapi	A 'Whoapi' API Client
WhopGenome	High-Speed Processing of VCF, FASTA and Alignment Data
widals	Weighting by Inverse Distance with Adaptive Least Squares for Massive Space-Time Data
widenet	Penalized Regression with Polynomial Basis Expansions
wikibooks	Functions and datasets of the german WikiBook "GNU R"
WikidataR	API Client Library for 'Wikidata'
WikipediaR	R-Based Wikipedia Client
wikipediatrend	Public Subject Attention via Wikipedia Page View Statistics
WikipediR	A MediaWiki API Wrapper
WikiSocio	A MediaWiki API Wrapper
WilcoxCV	Wilcoxon-based variable selection in cross-validation
wildlifeDI	Calculate Indices of Dynamic Interaction for Wildlife Telemetry Data
wildpoker	Best Hand Analysis for Poker Variants Including Wildcards
windex	windex: Analysing convergent evolution using the Wheatsheaf index
wingui	Advanced Windows Functions
wiod	World Input Output Database 1995-2011
WiSEBoot	Wild Scale-Enhanced Bootstrap
withr	Run Code 'With' Temporarily Modified Global State
wkb	Convert Between Spatial Objects and Well-Known Binary Geometry
wle	Weighted Likelihood Estimation
WMCapacity	GUI Implementing Bayesian Working Memory Models
WMDB	Discriminant Analysis Methods by Weight Mahalanobis Distance and bayes
wmlf	Wavelet Leaders in Multifractal Analysis
wmtsa	Wavelet Methods for Time Series Analysis
wnominate	WNOMINATE Roll Call Analysis Software
woe	Computes Weight of Evidence and Information Values
wordbankr	Accessing the Wordbank Database
wordcloud	Word Clouds
wordmatch	Matches words in one file with words in another file
wordnet	WordNet Interface
WordPools	Classical word pools used in studies of learning and memory
wPerm	Permutation Tests
wpp2008	World Population Prospects 2008
wpp2010	World Population Prospects 2010
wpp2012	World Population Prospects 2012
wpp2015	World Population Prospects 2015
wppExplorer	Explorer of World Population Prospects
wq	Exploring Water Quality Monitoring Data
wqs	Weighted Quantile Sum Regression
wrassp	Interface to the ASSP Library
WrightMap	IRT Item-Person Map with 'ConQuest' Integration
write.snns	Function for exporting data to SNNS pattern files
WriteXLS	Cross-Platform Perl Based R Function to Create Excel 2003 (XLS) and Excel 2007 (XLSX) Files
WRS2	A Collection of Robust Statistical Methods
wrspathrow	Functions for working with Worldwide Reference System (WRS)
wrspathrowData	Data used by the wrspathrow package
wrswoR	Weighted Random Sampling without Replacement
wrswoR.benchmark	Benchmark and Correctness Data for Weighted Random Sampling Without Replacement
wru	Who Are You? Bayesian Prediction of Racial Category Using Surname and Geolocation
wskm	Weighted k-Means Clustering
wsrf	Weighted Subspace Random Forest for Classification
wSVM	Weighted SVM with boosting algorithm for improving accuracy
wtcrsk	Competing Risks Regression Models with Covariate-Adjusted Censoring Weights
WufooR	R Wrapper for the 'Wufoo.com' - The Form Building Service
wux	Wegener Center Climate Uncertainty Explorer
WWGbook	Functions and datasets for WWGbook
0	
x12	x12 - wrapper function and structure for batch processing
x12GUI	X12 - Graphical User Interface
x13binary	Provide the 'x13ashtml' Seasonal Adjustment Binary
XBRL	Extraction of Business Financial Information from 'XBRL' Documents
x.ent	eXtraction of ENTity
xergm	Extensions of Exponential Random Graph Models
xergm.common	Common Infrastructure for Extensions of Exponential Random Graph Models
xgboost	Extreme Gradient Boosting
xgobi	Interface to the XGobi and XGvis programs for graphical data analysis
xhmmScripts	XHMM R scripts
XHWE	X Chromosome Hardy-Weinberg Equilibrium
XiMpLe	A Simple XML Tree Parser and Generator
xkcd	Plotting ggplot2 Graphics in an XKCD Style
XLConnect	Excel Connector for R
XLConnectJars	JAR dependencies for the XLConnect package
xlsx	Read, write, format Excel 2007 and Excel 97/2000/XP/2003 files
xlsxjars	Package required POI jars for the xlsx package
xmeta	A Toolbox for Multivariate Meta-Analysis
Xmisc	Xiaobei's miscellaneous classes and functions
XML	Tools for Parsing and Generating XML Within R and S-Plus
xml2	Parse XML
XML2R	EasieR XML data collection
XMRF	Markov Random Fields for High-Throughput Genetics Data
XNomial	Exact Goodness-of-Fit Test for Multinomial Data with Fixed Probabilities
xoi	Tools for Analyzing Crossover Interference
xpose4	Tools for Nonlinear Mixed-Effect Model Building and Diagnostics
xseq	Assessing Functional Impact on Gene Expression of Mutations in Cancer
xtable	Export Tables to LaTeX or HTML
xtal	Crystallization Toolset
xtermStyle	Terminal Text Formatting Using Escape Sequences
xts	eXtensible Time Series
xVA	Calculates Credit Risk Valuation Adjustments
xyloplot	A Method for Creating Xylophone-Like Frequency Density Plots
yacca	Yet Another Canonical Correlation Analysis Package
yaImpute	Nearest Neighbor Observation Imputation and Evaluation Tools
yakmoR	A Simple Wrapper for the k-Means Library Yakmo
YaleToolkit	Data exploration tools from Yale University
yaml	Methods to convert R data to YAML and back
ycinterextra	Yield curve or zero-coupon prices interpolation and extrapolation
yCrypticRNAs	Cryptic Transcription Analysis in Yeast
yhat	Interpreting Regression Effects
yhatr	R Binder for the Yhat API
YieldCurve	Modelling and estimation of the yield curve
ykmeans	K-means using a target variable
YplantQMC	Plant Architectural Analysis with Yplant and QuasiMC
YPmodel	The Short-Term and Long-Term Hazard Ratio Model for Survival Data
YuGene	A Simple Approach to Scale Gene Expression Data Derived from Different Platforms for Integrated Analyses
yuima	The YUIMA Project Package for SDEs
yummlyr	R Bindings for Yummly API
zCompositions	Imputation of Zeros and Nondetects in Compositional Data Sets
ZeBook	ZeBook Working with dynamic models for agriculture and environment
Zelig	Everyone's Statistical Software
ZeligChoice	Zelig Choice Models
zendeskR	Zendesk API Wrapper
zetadiv	Functions to Compute Compositional Turnover Using Zeta Diversity
zic	Bayesian Inference for Zero-Inflated Count Models
ZillowR	R Interface to Zillow Real Estate and Mortgage Data API
ZIM	Zero-Inflated Models for Count Time Series with Excess Zeros
zipcode	U.S. ZIP Code database for geocoding
zipfR	Statistical models for word frequency distributions
zoeppritz	Seismic Reflection and Scattering Coefficients
zoib	Bayesian Inference for Beta Regression and Zero-or-One Inflated Beta Regression
zoo	S3 Infrastructure for Regular and Irregular Time Series (Z's Ordered Observations)
zooaRch	Analytical Tools for Zooarchaeological Data
zooimage	Analysis of numerical zooplankton images
zoom	A spatial data visualization tool
zoon	Reproducible, Accessible & Shareable Species Distribution Modelling
ZRA	Dynamic Plots for Time Series Forecasting
ztable	Zebra-Striped Tables in LaTeX and HTML Formats
zyp	Zhang + Yue-Pilon trends package

This data is recorded every 10 seconds in real-time. 
 There are two human-readable timestamps, and a nanoseconds-since-epoch timestamp. The latter is best to use, and it should match the second human-readable timestamp. The time is set to UTC timezone. the r code should do the following Tasks: 
 
 1. Which instrument (any column with *_PNL" it int's header) has the highest instrument P/L on the day, and what is the value?
 2. Plot the relationship between instrument P/L and instrument volumes using the final values for the day 
 3. there are some erroneous data in the timestamps of the log file. Find and fix the data  --  no issues found
 4. Plot the portfolio P&L by timestamp (so X-axis is nanoseconds since epoch, Y-axis is portfolio cumulative P/L)
 
  date <- read.csv("C:\\Users\\devak\\Desktop\\Freelancer\\samplebig\\big.csv")

1. Added columns HumanTimestamp(D), Match(E)
2. D==CONCATENATE(TRIM(B2)," ",TRIM(C2))
3. E=EXACT(A2,D2)
4. Choose all records with E=1

bg = read.csv(file="big.csv",head=TRUE,sep=",")  

	myday       A_PNL    A_PHL    B_PNL    B_PHL    C_PNL
  2016-02-01    1        3        2        3        3
  2016-02-01    2        2        1        2        1
  2016-02-01    0        3        0        3        4
  2016-02-02    1        3        2        3        1
  2016-02-02    2        3        1        3        0
  2016-02-02    1        3        1        3        1
  2016-02-03    1        2        2        2        0
  2016-02-03    2        3        1        3        0
  2016-02-03    1        3        1        3        0

Aggregrate

colnames(X)[c(1,2)] <- c("good", "better")  -- change 1,2 column name
colnames(test)[c(1)] <- c("manoj")  -- change 1 column
colnames(test) <- c("a","b","c") -- change all three column


id <- c(1,1,1,2,2,2,3)
sal<- c(2,3,4,1,5,3,6)
test <- data.frame(id,sal)
ag <- aggregate(sal~id, data=test, max)
ag <-aggregate(test, by=list(id), FUN=max, na.rm=TRUE)

merge(test,ag)

ag2 <- aggregate(time.ms~subject, data=d2, max)

library(data.table)
test <- data.table(test)
test[, maxvalue:=max(sal), by = list(id)]


1	2
1	3
1       4
2	1
2	5
2	3
3	6

rbind - combine two data sets
  
In statistical modeling, regression analysis is a statistical process for estimating the relationships among variables

1) supervised learning - based on training data (input data). repeat until reached specified level of accuracy
eg:- spam detection, stock price at a time etc

	a) Logistic Regression -- stochastic (binary) - yes/no - Wald test - p-value 
	b) Back Propagation Neural Network.

2) Un-supervised learning - 
eg:-  A model is prepared by deducing structures present in the input data. This may be to extract general rules
	.It may through a mathematical process to systematically reduce redundancy, or it may be to organize data by similarity
	
	a) Apriori algorithm -- Applied in Transactional Databases
	b) K-Means -- 

3) 	
Statistics/Mathematics

Pearson Coefficient - r 

r(x,y) = COV(X,Y)/( Stdev(X) * Stdev(Y) )
r = 1 , +ve correlation
r = -1, -ve correlation
r = 0, no correrlation

1
Contents
Revision History ................................................................................................................. 3
Rules by Country ................................................................................................................ 3
1. WBCMORST: Originator Structuring ..................................................................... 5
2. WBCMBEST: Beneficiary Structuring ................................................................... 6
3. WBCMOCOB: Originator Country Different From Ordering Bank Country ......... 7
4. WBCMBCBB: Beneficiary Country Different From Beneficiary Bank Country ... 8
5. WBBMRDAM: Beneficiary with Multiple Round Dollar Amounts ....................... 9
6. WBOMRDAM: Originator with Multiple Round Dollar Amounts ....................... 10
7. WBCMOSDT: Small Dollar Transactions ............................................................ 11
8. GRROUNDT: Roundtrips Monthly Aggregation .................................................. 12
9. GRTRALPG: Transactions Above Limit (Amount) Peer Group Analysis ........... 13
10. HKCSBROS: Boiler Room Scam .......................................................................... 14
11. GRLREACC: Lately Reactivated Accounts .......................................................... 15
12. GRHRCJIF: Transactions From OR To High Risk Countries / Jurisdictions -
Internal Flag ..................................................................................................................... 16
13. WBCMOMCJ: The routing of Transactions through Multiple Countries /
Jurisdictions ...................................................................................................................... 17
14. WBCMSOSB/PBCMSOSB: Multiple Transactions Same Originator and Same
Beneficiary ....................................................................................................................... 18
15. GRRNDSUM: Round Sum Transactions .............................................................. 20
16. GROFFSHR - To and From Off-Shore Company ................................................. 21
17. WBCMMOSB: Multiple Transactions - Multiple Originators and Same
Beneficiary ....................................................................................................................... 23
18. WBCMSOMB: Multiple Transactions -Same Originator and Multiple
Beneficiaries ..................................................................................................................... 24
19. GRMTLIWL: Monthly Internal Watch List .......................................................... 25
20. PBAGCATX: Aggregate Large Cash Transactions ............................................... 27
21. GRTRVRCU: Transactions Beyond x% Variance ................................................ 28
22. CBCHLFDW: Large and Frequent Deposits / Withdrawals in a Charity .............. 29
23. SGNEWCUS: New Customer Monitoring for New To Bank Customers ............. 30
24. GRCASSTR: Cash Structuring .................................................................................. 31
25. GRCSDPWR: Deposit followed by numerous Wire transfers .............................. 32
26. GRMTLCAS: Monthly Cash Transaction Monitoring .......................................... 33
27. HKCASOCC: Cash Transactions for Profession codes ......................................... 34
28. HKCASRES: Cash Transactions for Residential codes ........................................ 35
29. GRCSSTBR: Structuring- Deposits at Various branches.......................................... 36
2
30. GRPPDTXN - PUPID ............................................................................................ 37
31. PBCASTXN: Cash transactions for Dynamic Risk review ........................................ 38
32. GRWKWASH: Rapid Movement of funds – Weekly Wash ................................. 39
33. GRLOTSCM: Lottery Scam .................................................................................. 40
34. GRHREIWL: High Risk Entity Match- SCB Internal Watch List. ....................... 41
35. WBCMOBHR: Multiple Transactions – Originator with Beneficiaries in High
Risk Countries/Jurisdictions ............................................................................................. 42
36. WBCMBOHR: Multiple Transactions – Beneficiary with Originators in High Risk
Countries/Jurisdictions ..................................................................................................... 43
37. WBCMBWSO: Multiple Transactions – Beneficiary with Same Originator ........ 44
38. WBCMOWSB: Multiple Transactions – Originator with Same Beneficiary ........ 45
39. WBCMOMTC: Originator with Multiple Transactions through Multiple Countries
46
40. WBCMBMTC: Beneficiary with Multiple Transactions through Multiple
Countries .......................................................................................................................... 48
41. WBBTOCTH: Beneficiary with Transactions above Country Thresholds ........... 50
42. WBOTOCTH: Originator with Transactions above Country Thresholds ............. 52
43. WBBSTOCT: Beneficiary with Single Transaction above Country Threshold .... 54
44. WBOSTOCT: Originator with Single Transaction above Country Threshold ...... 55
45. WBOHRBHR: Originator in High Risk Country/Jurisdiction with Transactions
above Country Thresholds with Beneficiaries in High Risk ............................................ 56
46. WBBHROHR: Beneficiary in High Risk Country/Jurisdiction with Transactions
above Country Thresholds from Originators in High Risk Countries/Jurisdictions ........ 58
47. WBCMRTOB: Roundtrips Between an Originator and Beneficiary ..................... 60
48. WBOST6MA: Originator Structuring – 6 Month Aggregation ............................. 61
49. WBBST6MA: Beneficiary Structuring – 6 Month Aggregation ........................... 63
50. WBSOSB6MA: Multiple Transactions - Same Originator/Same Beneficiary – 6
Month Aggregation .......................................................................................................... 65
51. WBOMT2MA: Originator with Multiple Transactions – 2 Month Aggregation .. 67
52. WBOMT6MA: Originator with Multiple Transactions – 6 Month Aggregation .. 68
53. WBBMT2MA: Beneficiary with Multiple Transactions – 2 Month Aggregation 70
54. WBBMT6MA: Beneficiary with Multiple Transactions – 6 Month Aggregation 72
55. SGTRVRCU: Transactions Beyond x% Variance ................................................. 73
56. USG01THC : Transactions Above Threshold – Monthly ..................................... 74
57. PBNEWCUS: New Customer Monitoring ............................................................ 76
58. GRWASHTR: Wash Transactions at a Customer Level ....................................... 77
59. CBCASTXN: Cash Transaction Monitoring ......................................................... 78
60. GLCASBLR: Boiler Room Scam for Cash ........................................................... 79
61. INIBACSA: Cash Deposit Structuring A .............................................................. 80
3
62. INIBACSB: Cash Deposit Structuring B ............................................................... 81
63. GRTRNEKY: New Account Monitoring for First 6 Months ................................ 82
64. WBCSPDMB: Cheques Same Payee Multiple Banks ........................................... 84
65. WBSQCKNO – Same payee Multiple sequential cheques .................................... 86
66. GRHVLTXN: Surge in Volume Transactions ....................................................... 88
67. UKLREACC: UK Lately Reactivated Accounts ................................................... 89
68. GRORGBLR: Boiler Room Scam - Originator ..................................................... 90
69. HKTALAMT – Transactions above Limit (Amount) ............................................ 91
Glossary of Terms ............................................................................................................ 92
Addenda China ................................................................................................................. 93
Addenda Germany ............................................................................................................ 94
Addenda Hong Kong ........................................................................................................ 95
Addenda India .................................................................................................................. 96
Addenda Indonesia ........................................................................................................... 97
Addenda Japan ................................................................................................................. 98
Addenda Jersey ................................................................................................................. 99
Addenda Malaysia .......................................................................................................... 100
Addenda Nigeria ............................................................................................................. 101
Addenda Pakistan ........................................................................................................... 102
Addenda Philippines ....................................................................................................... 103
Addenda Singapore ........................................................................................................ 104
Addenda Switzerland ..................................................................................................... 105
Addenda Taiwan ............................................................................................................. 106
Addenda Thailand .......................................................................................................... 107
Addenda UAE ................................................................................................................ 108
Addenda UK - PvB ......................................................................................................... 109
Addenda UK - WB ......................................................................................................... 110
Addenda US ................................................................................................................... 111
Revision History
Versio
n
Author Date Change Description
1.0 SCB •
1.1 Paresh Uppal 30-09-2014 • Adding the missing rules and their pseudo
codes
Rules by Country
4
S.NO DS rules
China
Germany
Hong Kong
India
Indonesia
Japan
Jersey
Malaysia
Nigeria
Pakistan
Philippines
Singapore
Switzerland
Taiwan
Thailand
UAE
UK - PvB
UK - WB
US
1 GRHREIWL High Risk Entity SCB Internal Watch List      
2 HKTALAMT Transactions Above Limit Amount  
3 UKLREACC UK Lately Reactivated Accounts 
4 GRPPDTXN PUPID Transactions                   
5 CBCASTXN Cash Transactions Monitoring  
6 GRCSSTBR Structuring Deposits at Various Branches              
7 HKCASOCC Cash Transactions Monitoring For Profession Codes           
8
HKCASRES Cash Transactions Monitoring For Residential Codes 
9 GRWKWASH Rapid Movement of funds Weekly Wash                 
10 PBAGCATX Large Cash Transactions Monitoring               
11 PBNEWCUS New Customer Monitoring      
12
GRTRVRCU
Transactions Beyond 10% Variance at Customer
Level
                
13 GRWASHTR Wash Transactions at Customer Level      
14
CBCHLFDW
Large and Frequent Deposits_Withdrawals in a
Charity
                
15
SGTRVRCU
Transactions Beyond 10% Variance at Customer
Level

16 SGNEWCUS Singapore New Customer Monitoring            
17 GRCASSTR Cash Structuring              
18 GRCSDPWR Deposits followed by Wire Transfers               
19 GRCASBLR Boiler Room Scam for Cash 
20 GRMTLCAS Monthly Cash Transaction Monitoring           
21 WBCMRTOB Round Trips between Originator and Beneficiary 
22
WBCMOBHR
Multiple Transactions Originator with Beneficiaries
in High Risk Countries_Jurisdictions
 
23
WBCMBOHR
Multiple Transactions Beneficiary with Originators in
High Risk Countries_Jurisdictions
 
24
WBCMSOSB
Multiple Transactions Same Originator Same
Beneficiary
   
25
WBCMBWSO
Multiple Transactions Beneficiary With Same
Originator
 
26 WBCMOWSB Originator With Same Beneficiary  
27 WBCMORST Originator Structuring                   
28 WBCMBEST Beneficiary Structuring                   
29
WBCMOCOB
Originator Country Different From Ordering Bank
Country
                  
30
WBCMOMTC
Originator with Multiple Transactions Through
Multiple Countries
 
31
WBCMBMTC
Beneficiary with Multiple Transactions Through
Multiple Countries
 
32
WBCMBCBB
Beneficiary Country Different From Beneficiary Bank
Country
                  
33 WBBMRDAM Beneficiary With Multiple Round Dollar Amounts          
34 WBOMRDAM Originator With Multiple Round Dollar Amounts          
35 WBBTOCTH Beneficiary Txns above Country Threshold  
36 WBOTOCTH Originator Txns above Country Threshold  
37 WBBSTOCT Beneficiary Single Transaction over Threshold  
38 WBOSTOCT Originator Single Transaction over Threshold  
39
WBOHRBHR Originator With High Risk to High Risk Over Threshold  
40
WBBHROHR
Beneficiary With High Risk to High Risk Over
Threshold
 
41 WBOST6MA Originator Structuring 6 Month Aggregation 
42 WBBST6MA Beneficiary Structuring 6 Month Aggregation 
43
WBSOSB6MA
Same Originator Same Beneficiary 6 Month
Aggregation

44
WBOMT2MA
Originator with Multiple Transactions 2 Month
Aggregation

45
WBOMT6MA
Originator with Multiple Transactions 6 Month
Aggregation

46
WBBMT2MA
Beneficiary with Multiple Transactions 2 Month
Aggregation

47
WBBMT6MA
Beneficiary with Multiple Transactions - 6 Month
Aggregation

48 WBCMOSDT Small Dollar Transactions         
49 GRROUNDT Round Trips Monthly Aggregation                  
50
GRTRALPG
Transactions Above Limit Amount Peer Group
Analysis
                  
51 HKCSBROS Boiler Room Scam                
52 GRTRNEKYC New Account Monitoring For First Six Months 
53 GRLREACC Lately Reactivated Accounts                  
54
GRHRCJIF
Transactions From_To High Risk
Countries_Jurisdictions Internal_Flag
                 
55 WBCMOMCJ Multiple Countries_Jurisdictions                  
56 GRHVLTXN Surge in Volume of Transactions 
57
PBCMSOSB
Multiple Transactions With the Same Originator and
Same Beneficiary
              
58 GRRNDSUM Round Sum Transactions                
59 GROFFSHR Off Shore Companies                   
60
WBCMMOSB
Multiple Transactions Multiple Originators and
Same Beneficiary
                
61
WBCMSOMB
Multiple Transactions Same Originator and Multiple
Beneficiaries
                
62 USG01THC Transactions above Threshold Monthly 
63 GRLOTSCM Lottery Scam 
64 GRMTLIWL Internal Watch list Monthly aggregation             
65 PBCASTXN Cash Transactions     
66 WBCSPDMB Cheques Same Payee Multiple Banks 
67 WBSQCKNO Cheques Same Payee Multiple Sequential Cheques 
68 GRORGBLR Boiler room Scam Originator 
69 INIBACSA Cash Deposit Structuring A 
70 INIBACSB Cash Deposit Structuring B 
30 31 33 32 28 22 31 28 27 27 27 28 33 27 27 30 29 21 36
5
1. WBCMORST: Originator Structuring
High Level Description
This scenario detects if the number of transactions for an Originator with an amount between
(Threshold1 & Threshold2) is equal to or greater than defined per month. When transactions for an
Originator meet these criteria, an Alert will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Sub Customer level.
Assumptions
The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
If the Originator field is blank, do not run the rule.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for an Originator is equal to or greater than defined per month
AND
The amount of each of the transactions is between Threshold1 and Threshold2
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Sub Customer Alert
Rule Code: WBCMORST
Rule Name: Originator Structuring
6
2. WBCMBEST: Beneficiary Structuring
High Level Description
This scenario detects if the number of transactions for a Beneficiary with an amount between
(Threshold1 & Threshold2) is equal to or greater than defined per month. When transactions for a
Beneficiary meet these criteria, an Alert will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Sub Customer level.
Assumptions
The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
If the Beneficiary field is blank, do not run the rule.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for a Beneficiary is equal to or greater than defined per month
AND
The amount of each f the transactions is between Threshold1 and Threshold2
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Sub Customer Alert
Rule Code: WBCMBEST
Rule Name: Beneficiary Structuring
7
3. WBCMOCOB: Originator Country Different From Ordering Bank
Country
High Level Definition
This scenario detects if the Originator‘s Country differs from the Ordering Bank‘s Country >= defined
transactions in a month with an aggregate amount equal to/greater than a Threshold1 amount
(excluding transactions with amounts less than Threshold2). When transactions for an Originator meet
these criteria, an Alert will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Sub Customer level.
Assumptions
• The Originator Country is identified by the following logic:
Transaction-daily_yyyymmdd.originator_country_code
• The Ordering Bank Country is identified by the following logic:
Transaction-daily_yyyymmdd.ordering_bank_country_code
If a Country field is blank, do not run the rule.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The Originator Country Code is different from the Ordering Bank Country Code
AND
The total number of such transactions for an Originator is >= defined per month
AND
The amount of such transactions is equal to/greater than Threshold2
AND
The total dollar amount of the total number of such transactions for an Originator in a month
excluding transactions with a dollar amount less than Threshold2 is >= Threshold1
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit =”CB/SA/PB‘
THEN
Create Sub Customer Alert
Rule Code: WBCMOCOB
Rule Name: Originator Country Different from Ordering Bank Country
8
4. WBCMBCBB: Beneficiary Country Different From Beneficiary Bank
Country
High Level Definition
This scenario detects if the Beneficiary‘s Country differs from the Beneficiary Bank‘s Country on equal
to/greater than defined transactions in a month with an aggregate dollar amount equal to or greater
than Threshold1 (excluding transactions with amounts less than Threshold2). When transactions for a
Beneficiary meet these criteria, an Alert will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Sub Customer level.
Assumptions
• The Beneficiary Country is identified by the following logic:
Transaction-daily_yyyymmdd.beneficiary_country_code
• The Beneficiary Bank Country is identified by the following logic:
Transaction-daily_yyyymmdd.beneficiary_bank_country_code
If a Country field is blank, do not run the rule.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The Beneficiary Country Code is different from the Beneficiary Bank Country Code
AND
The total number of such transactions for a Beneficiary is >= defined per month
AND
The amount of such transactions is equal to/greater than Threshold2
AND
The total dollar amount of the total number of such transactions for a Beneficiary in a month
excluding transactions with a dollar amount less than Threshold2 is >= Threshold1
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit =”CB/SA/PB‘
THEN
Create Sub Customer Alert
Rule Code: WBCMBCBB
Rule Name: Beneficiary Country Different from Beneficiary Bank Country
9
5. WBBMRDAM: Beneficiary with Multiple Round Dollar Amounts
High Level Definition
This scenario detects if the number of transactions for a Beneficiary with a round amount equal
to/greater than defined threshold is equal to/greater than defined number of transactions per month.
When transactions for a Beneficiary meet these criteria, an Alert will be raised.
This rule is specific to Private Banking.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Sub Customer level.
Assumptions
• Round amount is defined as an amount ending in 3 zeros followed by 2 decimal place zeros.
E.g., n.000.00 or 3,000.00 but not 3,100.00 or 3,000.25
• The Beneficiary is identified by the following logic:
Transaction-daily_yyyymmdd.beneficiary
If the Beneficiary field is blank, do not run the rule.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The Beneficiary has transactions with a round amount equal to/greater than threshold
AND
The number of transactions is equal to/greater than defined number of transactions per month
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = “PB”
THEN
Create Sub Customer Alert
Rule Code: WBBMRDAM
Rule Name: Beneficiary with Multiple Round Dollar Amounts
10
6. WBOMRDAM: Originator with Multiple Round Dollar Amounts
High Level Definition
This scenario detects if the number of transactions for an Originator with a round amount equal
to/greater than defined threshold, is equal to/greater than defined number of transactions per month.
When transactions for an Originator meet these criteria, an Alert will be raised.
This rule is specific to Private Banking.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Sub Customer level.
Assumptions
• Round amount is defined as an amount ending in 3 zeros followed by 2 decimal place zeros,
E.g. n,000.00 or 3,000.00 but not 3,100.00 or 3,000.25
• The Originator is identified by the following logic:
Transaction_daily_yyyymmdd.originator
If the Originator field is blank, do not run the rule.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The Originator has transactions with a round dollar amount equal to/greater than defined
threshold
AND
The number of such transactions is equal to/greater than defined number of transactions per
month
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = “PB”
THEN
Create Sub Customer Alert
Rule Code: WBOMRDAM
Rule Name: Originator with Multiple Round Dollar Amounts
11
7. WBCMOSDT: Small Dollar Transactions
High Level Definition
This scenario detects if the number of electronic transfers for a Customer type for a specified product
with an amount between min threshold amount and maximum threshold is equal to/greater than
defined per month. When transactions for a Customer meet these criteria, an Alert will be raised.
This rule is specific to Private Banking.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Account level.
Assumptions
All Electronic transfers are available for interrogation mapped to TT AML Code – TT0009.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of Electronic transfers for an Account is equal to/greater than defined per month;
AND
The amount of each of the transactions is greater than min threshold and is less than max
threshold;
AND
Organization Unit = “PB”
THEN
Create Account Alert:
Rule Code:
Rule Name: Corporate - Small Dollar Transactions
12
8. GRROUNDT: Roundtrips Monthly Aggregation
High Level Description
The scenario is run monthly to detect when the transactions with an aggregate amount from an
Originator to a Beneficiary are matched to the aggregate amount of transactions from that Beneficiary
back to the Originator within +/-10% in the same month. The alert will have a code of ”GRROUNDT”
and a name of ”Round Trips-Monthly Aggregation”.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at account level.
Assumptions
• The originator is identified by the following logic:
Transaction_daily_yyyymmdd.originator_account_num or
Transaction_daily_yyyymmdd.debit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.originator
• The beneficiary is identified by the following logic:
Transaction_daily_yyyymmdd.beneficiary_account_num or
transaction_daily_yyyymmdd.credit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.beneficiary
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Transaction Aggregate amount 1 >= Threshold
AND
Transaction Aggregate Amount 1 is within +/-10% Transaction Aggregate Amount 2
AND
1 originator = 2 beneficiary
AND
1 beneficiary = 2 originator
AND
1 originator is NOT 1 beneficiary
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Account Alert
Rule Code GRROUNDT,
Rule Name Round Trips - Monthly Aggregation,
13
9. GRTRALPG: Transactions above Limit (Amount) Peer Group Analysis
High Level Definition
This scenario detects if actual monthly transaction amount is greater than/equal to the average
monthly amount based on peer group analysis.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Account Level
Assumptions
Peer group analysis based on ISIC codes and codes to be given by countries – the codes can be
added or deleted by the DS Manager in the front end.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Monthly txn amount aggregate >= Factor* avg monthly amount based on peer group
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Account Alert
Rule Code: GRTRALPG
Rule Name: Transactions Above Limit (Amount) Peer Group Analysis
14
10. HKCSBROS: Boiler Room Scam
High Level Definition
This scenario detects accounts whose country of incorporation is classified as high risk, the
incorporation date is less than 12 months, the count and amount of monthly incoming TT transactions
Is >= defined and Threshold1 and amount of monthly debit transactions is >= Threshold2. An alert is
generated in this case.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Account Level
Assumptions
List of Country of Incorporation that is High Risk would be provided.
ITT transactions are also mapped to TT0019 AML Code
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Country of Incorporation is High Risk AND ISIC Code is within defined
AND
Account Since Opened Date < 12 months
AND
Count of monthly Credit TT
2
Transactions is >= defined
AND
Amount of monthly Credit TT transactions >= Threshold1
AND
Amount of monthly Debit Transactions >= Threshold2
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Account Alert
Rule Code: HKCSBROS
Rule Name: Boiler Room Scam
15
11. GRLREACC: Lately Reactivated Accounts
High Level Definition
This scenario detects reactivated accounts that have had more than defined transactions in the last
month from the reactivation date.
Execution Frequency
This detection scenario will be executed on a Daily basis.
Aggregation Level
The data will be monitored at the Account Level
Assumptions
It is assumed that the reactivation date would be available which could be then used to compare if the
reactivation happened in the last 60 days.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Account Reactivation Date < 60 days (soft coded)
AND
Number of Transactions > Transaction count
AND
Total amount of transactions > Threshold amt
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Account Alert
Rule Code: GRLREACC
Rule Name: Lately Reactivated Accounts
16
12. GRHRCJIF: Transactions From OR To High Risk Countries /
Jurisdictions - Internal Flag
High Level Description
This scenario detects when a payment transaction is from a country with an internal risk flag of 'Y' OR
to a country with an internal risk flag of 'Y'.
This scenario includes a first level segmentation for customer versus bank to bank transactions.
When the transaction is to or from a country with an internal risk flag of ”Y‘ , and the country codes are
not the same (e.g. intra-country transfers within one country or jurisdiction) an alert will be raised in the
case management module. The alert will have a code of ”GRHRCJIF” and a name of ”Transactions
from or to High Risk Countries/Jurisdictions -Internal Flag‘.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at transaction level.
Assumptions
Standard Chartered Bank will populate the internal high-risk flag on the country table to indicate which countries
are high risks.
From will be defined as the following country codes:
Transaction_daily_yyyymmdd.originator_country_code (if not a HR country check next)
Transaction_daily_yyyymmdd.ordering_bank_country_code (if not a HR country check next)
Transaction_daily_yyyymmdd.sending_bank_country_code (if not a HR country check next)
Transaction_daily_yyyymmdd.debit_party_country_code (if not a HR country check next)
Transaction_daily_yyyymmdd.Correspondent bank_country_code (if not a HR country check next)
Transaction_daily_yyyymmdd.Intermediary bank_country_code
To will be defined as the following country codes:
Transaction_daily_yyyymmdd.beneficiary_country_code (if not a HR country check next)
Transaction_daily_yyyymmdd.beneficiary_bank_country_code (if not a HR country check next)
Transaction_daily_yyyymmdd.credit_party_country_code
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Transaction is to OR from a High Risk Country
AND
Transaction is a payment message OR WIRE Transaction (103)
AND
Bank to Bank Flag = N or Blank
AND
Transaction Amount > Threshold AND Organization Unit = ‘WB/CB/SA/PB‘
THEN
Aggregate transactions
AND IF
No of such aggregated transactions>=Transaction Count Threshold
THEN
Create Alert
Rule Code: GRHRCJIF,
Rule Name: Transactions from or to High Risk Countries/Jurisdictions -Internal Flag,
17
13. WBCMOMCJ: The routing of Transactions through Multiple
Countries / Jurisdictions
High Level Description
This scenario detects the routing of Payment transactions (MT103) through several countries or
jurisdictions.
This scenario includes a first level segmentation for customer versus bank to bank transactions.
When the transactions are routed through several jurisdictions an alert will be raised in the case
management software. The alert will have a code ”WBCMOMCJ” and a name of ”Multiple
Countries/Jurisdictions‘.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at transaction level.
Assumptions
• Multiple jurisdictions are identified by examining the originator_country,
ordering_bank_country, sending_bank_country, debit_party_country, credit_party_country,
intermediary_bank_country, beneficiary_bank_country, beneficiary_country,
correspondent_bank_country codes of a customer’s customer‘s transaction
• A country field which is Blank or Null is not to be counted as part of the threshold >= 4
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Count of DISTINCT Country Code in (originator_country, ordering_bank_country,
sending_bank_country, debit_party_country, credit_party_country, intermediary_bank_country,
beneficiary_bank_country, beneficiary_country, correspondent_bank_country) is equal to or greater
than 4
AND
Bank to Bank Flag = N or Blank
AND
(Transaction is a Payment transaction (MT103)
OR
Transaction type is WIRE)
AND
Transaction Amount >= Threshold
THEN
Aggregate transactions
AND IF
No. of such aggregated transactions >= Transaction count threshold
AND Organization Unit =”WB/CB/SA/PB‘
THEN
Create Account Alert
Rule Code WBCMOMCJ,
Rule Name Multiple Countries/Jurisdictions,
18
14. WBCMSOSB/PBCMSOSB: Multiple Transactions Same Originator
and Same Beneficiary
High Level Description
This scenario detects the multiple transactions with the same originator and the same beneficiary.
This scenario includes a first level segmentation for customer versus bank to bank transactions.
When multiple transactions have the same originator and the same beneficiary an alert will be raised
in the case management software. The alert will have a code of ’WBCMSOSB’/’PBCMSOSB‘ and a
name of ”Multiple Transactions oe Same Originator/Same Beneficiary‘.
Execution Frequency
This detection scenario will be executed on a monthly basis
Aggregation Level
The data will be monitored at Sub Customer level for Incoming wires and Account level for Outgoing
wires.
Assumptions
US- KYCC
• The originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
Revised KYCC – Other than US- KYCC
The originator is identified by the following logic:
Transaction_daily_yyyymmdd.originator_account_num or
Transaction_daily_yyyymmdd.debit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.originator
The beneficiary is identified by the following logic:
Transaction_daily_yyyymmdd.beneficiary_account_num or
transaction_daily_yyyymmdd.credit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.beneficiary
WBCMSOSB only works with the TMW US KYCC package as it is looking for transactions where the
originator is the same as the beneficiary i.e. the same entity. The Revised KYCC package does not
allow for the originator to be the same as beneficiary. WBCMSOSB fires always on the Originator and
is the logic to be used for clearing centers.
A new rule PBCMSOSB has been written to cater for the revised KYCC package where the DS is
looking for the same entity i.e. the originator is the same as the beneficiary but it is looking for a set of
transactions where there are a number of transactions having the same originator and the same
beneficiary.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Number of same originator and same beneficiary transactions is >= defined
AND
Bank to Bank Flag = N or Blank
AND
Transaction Amount aggregate is >= Threshold
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Sub Customer/Account Alert
19
Rule Code WBCMSOSB, PBCMSOSB
Rule Name Multiple Transactions - Same Originator/Same Beneficiary
20
15. GRRNDSUM: Round Sum Transactions
High Level Definition
This scenario detects if the number of round sum transactions in an account is more than defined
count.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Account level.
Assumptions
• Round Sum amount is defined as an amount ending in zeros followed by 2 decimal place
zeros, e.g. n,000.00
• This rule will monitor all customer transactions unless specified.
• Bank to Bank transactions and Book transfers are not part of interrogation by this rule.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The transaction amount ends with number of zeroes >= specified count
AND
The total number of such transactions > specified transaction count
AND
Bank to Bank Flag = Y, N or Blank;
AND
Organization Unit = “CB/WB/SA”
THEN
Create Account Alert:
Rule Code:
Rule Name: Round Sum Transactions
21
16. GROFFSHR - To and From Off-Shore Company
High Level Description
Alerts to be generated whenever the offshore country table (List provided by in-country) matches and
the transaction is above the threshold and the number of transaction is also above threshold.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Account level.
Assumptions
• The DS is expected to monitor only Cross border transactions identified by TT0019, TT0020.
• Any transaction types that should be excluded must be populated in the agg_parameters table
with the type set to ‘TXN’.
• Only transactions with the bank_to_bank_flag set to ‘N’ will be included.
• Only transactions with the curr_ref_calc_flag set to ‘Y’ will be included.
• Country_risk_flags table need to be updated and set the column offshore_fin_center_flag = ‘Y’
for all applicable Off shore countries
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions
WHERE
Bank_to_bank_flag = ‘N’
AND
Curr_ref_calc_flag = ‘Y’
AND
(Beneficiary Country Code,
OR
Originator Country Code,
OR
Beneficiary Bank Country Code,
OR
Ordering Bank Country Code,
OR
Sending Bank Country Code,
OR
Debit Party Country Code,
OR
Credit Party Country Code,
OR
Correspondent Bank Country Code,
OR
Intermediary Bank Country Code)
Matches an entry in the country_risk_flags
WHERE
offshore_fin_center_flag = ‘Y’
AND
Transaction amount > XXXXXX
22
Is greater than the defined transaction count
THEN
Create Account Alert
Rule Code: GROFFSHR
Rule Name: Transactions From/To Off Shore countries
23
17. WBCMMOSB: Multiple Transactions - Multiple Originators and
Same Beneficiary
High Level Description
This scenario detects the multiple transactions with different originators and the same beneficiary.
This scenario includes a first level segmentation for customer versus bank to bank transactions.
When multiple transactions have different originators identified by the count and the same beneficiary
an alert will be raised in the case management software. The alert will have a code of ”WBCMMOSB‘
and a name of ”Multiple Transactions -Multiple Originators and Same Beneficiary.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
The originator is identified by the following logic:
Transaction_daily_yyyymmdd.originator
The beneficiary is identified by the following logic
Transaction_daily_yyyymmdd.beneficiary
Transaction_daily_yyyymmdd.originator_account_num or
Transaction_daily_yyyymmdd.debit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.originator
The beneficiary is identified by the following logic:
Transaction_daily_yyyymmdd.beneficiary_account_num or
Transaction_daily_yyyymmdd.credit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.beneficiary
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Number of different originator and same beneficiary transactions in the month is > defined (No
of different Originators are > defined)
AND
Bank To Bank Flag = N or Blank
AND
Aggregate Base TXN Amount of the transactions is >= Threshold
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Sub Customer Alert
Rule Code WBCMMOSB,
Rule Name Multiple Transactions -Multiple Originators and Same Beneficiary,
24
18. WBCMSOMB: Multiple Transactions -Same Originator and Multiple
Beneficiaries
High Level Description
This scenario detects the multiple transactions with the same originator and different beneficiaries.
This scenario includes a first level segmentation for customer versus bank to bank transactions.
When multiple transactions have the same originator and different beneficiaries identified by the count
an alert will be raised in the case management software. The alert will have a code of ”WBCMSOMB‘
and a name of ”Multiple Transactions -Same Originator and Multiple Beneficiaries.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
The originator is identified by the following logic:
Transaction_daily_yyyymmdd.originator
The beneficiary is identified by the following logic:
Transaction_daily_yyyymmdd.beneficiary
Transaction_daily_yyyymmdd.originator_account_num or
Transaction_daily_yyyymmdd.debit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.originator
The beneficiary is identified by the following logic:
Transaction_daily_yyyymmdd.beneficiary_account_num or
transaction_daily_yyyymmdd.credit_party_account_num (if not exists, check next)
Transaction_daily_yyyymmdd.beneficiary
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Number of same originator and different beneficiary transactions in the month is > defined
(The number of different beneficiaries > defined)
AND
Bank To Bank Flag = N or Blank
AND
Aggregate Base TXN Amount of the transactions is >= Threshold
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Sub Customer Alert
Rule Code WBCMSOMB,
Rule Name Multiple Transactions -Same Originator and Multiple Beneficiaries,
25
19. GRMTLIWL: Monthly Internal Watch List
High Level Description
This scenario detects accounts which have transactions to OR from entities internally flagged as highrisk
entity.
When the transaction is to OR from an entity on an internal watch list an alert will be raised in the case
management software as a monthly alert for all transactions for that month in the account.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Account level.
Assumptions
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF Country has Unique Identifier logic defined then if the customer id of the transaction matches the
customer id in the Watch list.
OR
IF Country does not have unique identifier logic implemented then
(
Transaction is to OR from an originator on Watch list
OR
Transaction is to OR from an ordering_bank on Watch list
OR
Transaction is to OR from a sending_bank on Watch list
OR
Transaction is to OR from a debit_party on Watch list
OR
Transaction is to OR from a credit_party on Watch list
OR
Transaction is to OR from an intermediary_bank on Watch list
OR
Transaction is to OR from a correspondent_bank on Watch list
OR
Transaction is to OR from a beneficiary_bank on Watch list
OR
Transaction is to OR from a beneficiary on Watch list
OR
Transaction is to OR from an originating bank on Watch list
OR
Transaction is to OR from the BBI reference on Watch list
OR
Transaction is to OR from Senders Ref on Watch list
OR
Transaction is to OR from Related Reference on Watch list
OR
Transaction is to OR from Remittance Reference on Watch list
OR
Transaction is to OR from originator to beneficiary Reference on Watch list
)
AND
Organization Unit ‘WB/CB/SA/PB‘
THEN
26
Create Account Alert
Rule Code GRMTLIWL,
Rule Name Monthly Internal Watch List,
27
20. PBAGCATX: Aggregate Large Cash Transactions
High Level Definition
This scenario detects customers whose cash deposits or withdrawals are greater than or equal to
defined amount in pseudo code below in a month.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at the Customer Level
Assumptions
It is assumed that all accounts of a customer will be aggregated for obtaining total amount of cash
deposits / withdrawals in a month.
All cash transactions are mapped to TT0001.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
AND
Total Amount of Cash Withdrawals in a month >= Threshold
OR
Total Amount of Cash Deposits in a month >= Threshold
AND
Organization Unit =”PB‘
THEN
Create Customer Alert
Rule Code: PBAGCATX
Rule Name: Aggregate Large Cash Transactions Monitoring
28
21. GRTRVRCU: Transactions Beyond x% Variance
High Level Description
This scenario detects the variance between the actual transactions of a customer and the customer‘s
transactional activity based on past three months.
When the actual monthly transaction amount for a customer across all his accounts is greater than or
equal to a defined value and it is greater than a defined average monthly amount based on past 3
months for the customer who has a relationship for 4 months, an Alert will be raised.
Execution Frequency
This detection scenario will be run on a monthly basis.
Aggregation Level
The data will be monitored at customer level.
Assumptions
Anticipated Volume is loaded into system for generation of alerts based on expected volume. If not,
only the second condition would be executed.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
(Customer Acquisition date > 120 days
AND
Actual Monthly Txn Amt >= Threshold value
AND
Actual Monthly Txn Amt is >= Factor * Avg monthly amt based on past 3 months
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Customer Alert
Rule Code GRTRVRCU,
Rule Name Transactions Beyond x% Variance at Customer Level
29
22. CBCHLFDW: Large and Frequent Deposits / Withdrawals in a
Charity
High Level Definition
This scenario detects customers who hold charity accounts and the number and amount of
withdrawals/deposits in the last 30 days are more than what is defined.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at the Customer Level
Assumptions
It is assumed that charity accounts / customers can be identified by ISIC codes– the codes can be
added or deleted by the DS Manager in the front end.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Total Number of Withdrawals in last 30 days >= 10
AND
Total Amount of Withdrawals in last 30 days >= Threshold
OR
Total Number of Deposits in last 30 days >= 10
AND
Total Amount of Deposits in last 30 days >= Threshold
AND
Account is a Charity Account
AND
Organization Unit = ”WB/CB ‘
THEN
Create Customer Alert
Rule Code: CBCHLFDW
Rule Name: Large and Frequent Deposits / Withdrawals in a Charity
30
23. SGNEWCUS: New Customer Monitoring for New To Bank
Customers
High Level Definition
This scenario monitors New to Bank customers within the last 6 months and an alert is generated
when monthly or debit turnover is greater than or equal to Threshold amount.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at the Customer Level
Assumptions
New to Bank Customer in the last 6 months.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Customer Acquisition date < 6 months
AND
Total Debit in a Month >= Threshold
OR
Total Credit in a Month >= Threshold
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Customer Alert
Rule Code: SGNEWCUS
Rule Name: New Customer Monitoring
31
24. GRCASSTR: Cash Structuring
High Level Description
This scenario detects when a customer has more than defined structured cash transactions in his/her
accounts.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Customer level.
Assumptions
Cash transaction types available in the aggregation parameter table as TXN type
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Amount of each cash deposit or cash withdrawal transaction should be between Min threshold
and Max threshold
AND
Aggregated amount of such transactions >Aggregated threshold and the number of such
transactions > defined count then
THEN
Create Alert
Rule Code: GRCASSTR,
Rule Name: Cash Structuring,
32
25. GRCSDPWR: Deposit followed by numerous Wire transfers
High Level Description
This scenario detects when a customer has huge cash deposit which is followed by most of the fund
being wired out.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Customer level.
Assumptions
• Cash transaction types available in the aggregation parameter table as CASH type
• OTT Wire codes available in aggregation parameter table as WIRE type
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Number of CASH deposit transaction > Deposit count
AND
Amount of such CASH deposit transactions > Deposit Threshold
AND
Number of outward remittances > OTT count
AND
Total amount of such outward remittances > OTT Threshold
OR
Total aggregate amount of Cash deposits > Deposit Threshold
AND
Total amount of such outward remittances should be +/- defined percentage of Aggregate
Cash Deposit
THEN
Create Alert
Rule Code: GRCSDPWR,
Rule Name: Deposit followed by numerous Wire transfers,
33
26. GRMTLCAS: Monthly Cash Transaction Monitoring
High Level Description
This scenario detects customers who have cash deposits or withdrawals for at least a day greater than
or equal to Threshold amount.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Customer level.
Assumptions
• Cash transaction types need to be populated in aggregation parameters table.
• Thresholds need to be updated in the Scenario Values table.
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
There is at least one instance in a month
WHERE
Total Amount of Cash Withdrawals in a day >= Threshold1
OR
Total Amount of Cash Deposits in a day >= Threshold2
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Customer Alert
Rule Code: GRMTLCAS
Rule Name: Monthly Cash Transactions Monitoring
34
27. HKCASOCC: Cash Transactions for Profession codes
High Level Definition
This scenario detects low salaried customers whose cash deposits or withdrawals are greater than or
equal to Threshold amount in a day.
Execution Frequency
This detection scenario will be executed on a Daily basis.
Aggregation Level
The data will be monitored at the Customer Level
Assumptions
• It is assumed that all accounts of a customer will be aggregated for obtaining total amount
of cash deposits / withdrawals in a day.
• Low Salaried customers can be identified by Professional/Occupational codes. These can be
added by DS Manager in the front end.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Total Amount of Cash Withdrawals in a day >= Threshold1
OR
Total Amount of Cash Deposits in a day >= Threshold2
AND
Customers specified by Profession codes
AND
Organization Unit = ”CB ‘
THEN
Create Customer Alert
Rule Code: HKCASOCC
Rule Name: Cash Transactions Monitoring
35
28. HKCASRES: Cash Transactions for Residential codes
High Level Definition
This scenario detects customers for specific nationality whose cash deposits or withdrawals is greater
than or equal to Threshold amount in a day.
This rule is specific to Hong Kong.
Execution Frequency
This detection scenario will be executed on a Daily basis.
Aggregation Level
The data will be monitored at the Customer Level
Assumptions
It is assumed that all accounts of a customer will be aggregated for obtaining total amount of cash
deposits /withdrawals in a day.
Nationality for the customers for monitoring is available and added to the DS Rule logic.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Total Amount of Cash Withdrawals in a day >= Threshold1
OR
Total Amount of Cash Deposits in a day >= Threshold2
AND
Customers for specific nationality
AND
Organization Unit = ”HK CB “
THEN
Create Customer Alert
Rule Code: HKCASOCC
Rule Name: Cash Transactions Monitoring
36
29. GRCSSTBR: Structuring- Deposits at Various branches
High Level Description
This scenario detects when a customer has multiple cash deposits in various branches on the same
day.
Execution Frequency
This detection scenario will be executed on a Daily basis.
Aggregation Level
The data will be monitored at Customer level.
Assumptions
• Cash transaction types available in the aggregation parameter table as TXN type
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Number of Cash deposit transaction is > transaction count
AND
Amount of such transactions is > threshold
AND
Number of distinct branches > branch count then
THEN
Create Alert
Rule Code: GRCSSTBR,
Rule Name: Structuring- Deposits at Various branches,
37
30. GRPPDTXN - PUPID
High Level Description
Transactions sent by or to non-customers, also known as “Payable Upon Proper Identification”
(PUPID)
Execution Frequency
This detection scenario will be executed on a Daily basis.
Aggregation Level
The data will be monitored at Sub Customer level (Beneficiary).
Assumptions
• Any transaction types that should be excluded must be populated in the agg_parameters table
with the type set to ‘TXN’.
• Only transactions with the bank_to_bank_flag set to ‘N’ will be included.
• Only transactions with the curr_ref_calc_flag set to ‘Y’ will be included.
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
PUPID or PAY UPON is in the content of the below fields
1. TRANS_REF_DESC
2. ORIGINATOR_ACCOUNT_NUMBER
3. ORIGINATOR
4. BENEFICIARY_ACCOUNT_NUMBER
5. BENEFICIARY
6. BENEFICIARY_ADDRESS_1
7. BENEFICIARY_ADDRESS_2
8. BENEFICIARY_ADDRESS_3
9. BENEFICIARY_ADDRESS_4
10. TRANS_REF_DESC_4
THEN
Create Account Alert
Rule Code: GRPPDTXN
Rule Name: PUPID Transactions
38
31. PBCASTXN: Cash transactions for Dynamic Risk review
High Level Description
This scenario detects accounts which have more than specified cash deposits/withdrawals for a period
of past 12 months.
This scenario is specific to Private Banking.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at Customer level.
Assumptions
• Cash transaction types available in the aggregation parameter table as TXN type
High Level Pseudo code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Aggregate cash credits for the past 12 months>= threshold
OR
Aggregate cash debits for the past 12 months >= threshold
AND
Organization Unit = ”PB‘
THEN
Create Account Alert
Rule Code: PBCASTXN
Rule Name: Cash transactions for Dynamic Risk review
39
32. GRWKWASH: Rapid Movement of funds – Weekly Wash
High Level Definition
This scenario detects transactions where the amount that has been credited in to the accounts of a
customer is removed in full (+/-x%) in the past calendar week.
Execution Frequency
This detection scenario will be executed on a weekly basis.
Aggregation Level
The data will be monitored at the Customer Level
Assumptions
The threshold values are considered to be credit amount thresholds over past calendar week.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Amount Credited in past calendar week >= Threshold
AND
Amount debited in past calendar week is within +/-7% of Amount credited in past calendar
week)
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Customer Alert
Rule Code: GRWKWASH
Rule Name: Rapid Movement of funds – Weekly Wash
40
33. GRLOTSCM: Lottery Scam
High Level Definition
This scenario detects transactions having a pattern of Lottery Scam.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Account Level
Assumptions
The Wire IN and Cash transactions are populated in agg_parameters.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Nationality or Country of domicile is Taiwan/China
AND
Customer Acquisition date within 6 months
AND
Inward Remittance Transactions (TT0019) > threshold
AND
Cash withdrawal (TT0001) > threshold
AND
Organization Unit = ”CB‘
THEN
Create Account Alert
Rule Code: GRLOTSCM
Rule Name: Lottery Scam
41
34. GRHREIWL: High Risk Entity Match- SCB Internal Watch List.
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank. This scenario detects
transactions to/from entities internally flagged as high-risk entity and subject to monitoring.
When the transaction is to/from an entity on the SCB Internal Watch List an alert will be raised in the
case management software. The alert will have a code of ‘GRHREIWL’ and a name of ‘High Risk
Entity Match- SCB Internal Watch List’. The priority displayed on the end User screen will be set to 10
to indicate that the alert is High priority.
Execution Frequency
This detection scenario will be executed on a daily basis.
Aggregation Level
The data will be monitored at transaction level.
Assumptions
• The SCB Internal Watch List will be created and maintained by the Group List Manager using the
Norkom List Management software.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF (
Transaction is to/from an originator on Watch list OR
Transaction is to/from an ordering_bank on Watch list OR
Transaction is to/from a sending_bank on Watch list OR
Transaction is to/from a debit_party on Watch list OR
Transaction is to/from a credit_party on Watch list OR
Transaction is to/from an intermediary_bank on Watch list OR
Transaction is to/from a correspondent_bank on Watch list OR
Transaction is to/from a beneficiary_bank on Watch list OR
Transaction is to/from a beneficiary on Watch list
) AND
Organisation Unit = ‘US WB’
THEN
Create Transaction Alert
Rule Code GRHREIWL
Rule Name High Risk Entity Match - SCB Internal Watch List
Customer Customer Key
Account Account Key
Transaction Transaction Key
Priority 10
Risk Weighting 8
42
35. WBCMOBHR: Multiple Transactions – Originator with
Beneficiaries in High Risk Countries/Jurisdictions
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rules G3-10-O and G4-10-O. This scenario detects the multiple transactions where
an Originator sends funds to Beneficiaries in High Risk countries or jurisdictions multiple times in a
month. When an Originator sends funds to Beneficiaries in a High Risk country/jurisdiction more than
10 times in a month with a transaction amount equal to/greater than USD 100,000.00 an Alert will be
raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Originator country is identified by the following logic:
Transaction_daily_yyyymmdd.originator_country_code
• The Beneficiary country is identified by the following logic:
Transaction_daily_yyyymmdd.beneficiary_country_code
• The country risk is identified by the following logic: Internal_risk_flag = Y
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514 representing a
Book Transfer then count only a single occurrence of the Transaction Reference Number (TRN No.)
for either the transaction threshold or the dollar amount threshold.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF The Originator has Beneficiaries with country codes flagged as Internal risk flag = y
AND
Bank to Bank Flag = N or Blank
AND
Transaction Amount is equal to /greater than USD 100,000.00
AND
The number of transactions is greater than 10
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code WBCMOBHR
Rule Name Multiple Transactions –Originator with Beneficiaries in High Risk Countries
Customer Sub Customer Key,
Priority 5
Risk Weighting 7
43
36. WBCMBOHR: Multiple Transactions – Beneficiary with
Originators in High Risk Countries/Jurisdictions
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rules G3-10-B and G4-10-B. This scenario detects the multiple transactions where
a Beneficiary receives funds from Originators in High Risk countries or jurisdictions multiple times in a
month. When a Beneficiary receives funds from Originators in a High Risk country/jurisdiction more
than 10 times in a month, with a transaction amount equal to/greater than USD 100,000.00 an Alert
will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Originator country is identified by the following logic:
Transaction_daily_yyyymmdd.originator_country_code
• The Beneficiary country is identified by the following logic:
Transaction_daily_yyyymmdd.beneficiary_country_code
• The country risk is identified by the following logic: Internal_risk_flag = Y
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF
The Beneficiary has Originators with country codes flagged as Internal risk flag = Y
AND
Bank to Bank Flag = N or Blank
AND
Transaction Amount is equal to/greater than USD 100,000.00
AND
The number of transactions is greater than 10
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code WBCMBOHR
Rule Name Multiple Transactions – Beneficiary with Originators in High Risk
Customer Customer Key
Priority 5
Risk Weighting 7
44
37. WBCMBWSO: Multiple Transactions – Beneficiary with Same
Originator
High Level Description
This scenario detects the multiple transactions where a Beneficiary receives funds from the same
Originator multiple times in a month. This scenario is run for the Standard Chartered Bank’s GB
Wholesale Bank, and corresponds to the Look Back Project Rules G3-2-B and G4-2-B.
When a Beneficiary receives funds from the same Originator more than 9 times in a month for a total
amount greater than GBP 200,000.00 an Alert will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF The beneficiary has the same originator more than 9 transactions
AND
The total amount of the transactions is greater than or equal to GBP 100,000.00
AND
Customers Customer Flag = ‘Y’
AND
Organisation Unit = ‘GB WB’
THEN
Create Sub Customer Alert
Rule Code WBCMBWSO,
Rule Name Multiple Transactions –Beneficiary with same Originator,
Customer Sub Customer Key,
Priority 5,
Risk Weighting 7
45
38. WBCMOWSB: Multiple Transactions – Originator with Same
Beneficiary
High Level Description
This scenario detects the multiple transactions where an Originator sends funds to the same
Beneficiary multiple times in a month. This scenario is run for the Standard Chartered Bank’s US
Wholesale Bank, and corresponds to the Look Back Project Rules G3-2-O and G4-2-O. This scenario
includes a first level segmentation for customer versus bank to bank transactions. When an Originator
sends funds to the same Beneficiary more than 9 times in a month for a total amount greater than
GBP 100,000.00 an Alert will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
Pseudo Code
IF
The originator has the same beneficiary > 9 transactions
AND
The total amount of the transactions is equal to/greater than GBP 100,000.00
AND
Customers Customer Flag = ‘Y’
AND
Organisation Unit = ‘GB WB’THEN
THEN
Create Sub Customer Alert
Rule Code WBCMOWSB,
Rule Name Multiple Transactions – Originator with same Beneficiary,
Customer Sub Customer Key,
Priority 5,
Risk Weighting 7
46
39. WBCMOMTC: Originator with Multiple Transactions through
Multiple Countries
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-7-O and G4-7-O. This scenario detects if an Originator has 3 or more
transactions in a month each with a USD amount equal to/greater than USD 100,000.00, and each of
the transactions passed through multiple countries/jurisdictions. When transactions for an Originator
meet these criteria, an Alert will be raised.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Multiple Countries/Jurisdictions are defined by the following logic:
Transaction-daily_yyyymmdd.originator_country
Transaction-daily_yyyymmdd.ordering_bank_country
Transaction-daily_yyyymmdd.sending_bank_country
Transaction-daily_yyyymmdd.debit_party_country
Transaction_daily_yyyymmdd.credit_party_country
Transaction_daily_yyyymmdd.intermediary_bank_country
Transaction_daily_yyyymmdd.correspondent_bank_country
Transaction_daily_yyyymmdd.beneficiary_bank_country
Transaction_daily_yyyymmdd.beneficiary_country
• A Country Code field which is Blank or Null is not to be counted as part of the threshold = 4
• If the Originator field is blank, do not run the rule.
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF
The Originator has 3 or more transactions in a month
AND
The amount of each of the transactions is equal to/greater than USD 100,000.00
AND
The count of distinct Country Codes (as defined above) on each of the transactions is equal
to/greater than 4
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBCMOMTC
47
Rule Name: Originator with Multiple Transactions through Multiple Countries
Customer Sub Customer Key
Priority: 9
Risk Weighting: 9
48
40. WBCMBMTC: Beneficiary with Multiple Transactions through
Multiple Countries
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-7-B and G4-7-B. This scenario detects if a Beneficiary has 3 or more
transactions in a month each with a USD amount equal to/greater than USD 100,000.00, and each of
the transactions passed through multiple countries/jurisdictions. When transactions for a Beneficiary
meet these criteria, an Alert will be raised..
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Multiple Countries/Jurisdictions are defined by the following logic:
Transaction-daily_yyyymmdd.originator_country
Transaction-daily_yyyymmdd.ordering_bank_country
Transaction-daily_yyyymmdd.sending_bank_country
Transaction-daily_yyyymmdd.debit_party_country
Transaction_daily_yyyymmdd.credit_party_country
Transaction_daily_yyyymmdd.intermediary_bank_country
Transaction_daily_yyyymmdd.correspondent_bank_country
Transaction_daily_yyyymmdd.beneficiary_bank_country
Transaction_daily_yyyymmdd.beneficiary_country
• A Country Code field which is Blank or Null is not to be counted as part of the threshold = 4
• If the Beneficiary field is blank, do not run the rule.
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF
The Beneficiary has 3 or more transactions in a month
AND
The amount of each of the transactions is equal to/greater than USD 100,000.00
AND
The count of distinct Country Codes (as defined above) on each of the transactions is equal
to/greater than 4
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBCMBMTC
49
Rule Name: Beneficiary with Multiple Transactions through Multiple Countries
Customer Sub Customer Key
Priority: 9
Risk Weighting: 9
50
41. WBBTOCTH: Beneficiary with Transactions above Country
Thresholds
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-12-B and G4-12-B. This scenario detects if the total number of
transactions for a Beneficiary is greater than the total number of transactions threshold set for the
Beneficiary’s Country, and the total dollar amount of the Beneficiary’s transactions is greater than the
total dollar amount threshold set for the Beneficiary’s Country in a month. When transactions for a
Beneficiary meet these criteria, an Alert will be raised. For example, in the Look Back Project a
transaction threshold was derived for the UAE as a total of 13 transactions per month and a total dollar
amount threshold of equal to/greater than USD 1,000,000.00. If a Beneficiary in the UAE receives
more than 13 transactions in a month with a total dollar amount equal to/greater than USD
1,000,000.00 the Alert would generate.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Beneficiary Country is identified by the following logic:
Transaction-daily_yyyymmdd.beneficiary_country_code
• The Beneficiary Volume Threshold for each Country Code is defined in Table
• The Beneficiary Dollar Amount Threshold for each Country Code is defined in Table
• If the Beneficiary field is blank, do not run the rule.
• If the Beneficiary field is populated but the Beneficiary Country field is blank, default to the
Country Code = NQ representing Null.
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF
The total number of transactions for a Beneficiary in a month is greater than the total number
of transactions threshold set for the Beneficiary Country Code as defined in Table
AND
The total dollar amount of the total number of transactions for a Beneficiary in a month is
greater than the total dollar amount threshold set for the Beneficiary’s Country Code as
defined in Table
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBBTOCTH
51
Rule Name: Beneficiary with Transactions above Country Thresholds
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
52
42. WBOTOCTH: Originator with Transactions above Country
Thresholds
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-12-O and G4-12-O. This scenario detects if the total number of
transactions for an Originator is greater than the total number of transactions threshold set for the
Originator’s Country, and the total dollar amount of the Originator’s transactions is greater than the
total dollar amount threshold set for the Originator’s Country in a month. When transactions for an
Originator meet these criteria, an Alert will be raised.
For example, in the Look Back Project a transaction threshold was derived for the UAE as a total of 17
transactions per month and a total dollar amount threshold of equal to/greater than USD 1,000,000.00.
If an Originator in the UAE had more than 17 transactions in a month with a total dollar amount equal
to/greater than USD 1,000,000.00 the Alert would generate.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Originator Country is identified by the following logic:
Transaction-daily_yyyymmdd.originator_country_code
• The Originator Volume Threshold for each Country Code is defined in Table
• The Originator Dollar Amount Threshold for each Country Code is defined in Table
• If the Originator field is blank, do not run the rule.
• If the Originator field is populated but the Originator Country field is blank, default to the
Country Code = NQ representing Null.
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF
The total number of transactions for an Originator in a month is greater than the total number
of transactions threshold set for the Originator Country Code as defined in Table
AND
The total dollar amount of the total number of transactions for an Originator in a month is
greater than the total dollar amount threshold set for the Originator’s Country Code as defined
in Table
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
53
Create Sub Customer Alert
Rule Code: WBOTOCTH
Rule Name: Originator with Transactions above Country Thresholds
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
54
43. WBBSTOCT: Beneficiary with Single Transaction above Country
Threshold
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-11-B and G4-11-B. This scenario detects if a Beneficiary receives a
transaction with a USD amount greater than the single dollar amount threshold set for the Beneficiary’s
Country. When transactions for a Beneficiary meet these criteria, an Alert will be raised.
For example, in the Look Back Project a transaction threshold was derived for the UAE as a single
dollar amount threshold of USD 1,700,000.00. If a Beneficiary in the UAE receives a single transaction
with a USD amount greater than USD 1,700,000.00 the Alert would generate.
This Rule is designed for those entities (i.e. Beneficiaries), who have only a single transaction in a
month, and therefore are unlikely to hit the volume thresholds embedded within other Rules. This Rule
will capture large dollar amount transactions, and subject them to a review.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Beneficiary Country is identified by the following logic:
Transaction-daily_yyyymmdd.beneficiary_country_code (If Blank, check next)
• The Beneficiary Single Transaction Dollar Amount Threshold for each Country Code is
defined in Table
• If the Beneficiary field is blank, do not run the rule.
• If the Beneficiary field is populated but the Beneficiary Country field is blank, default to the
Country Code = NQ representing Null.
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF
A Beneficiary has a transaction in a month with a USD amount greater than the single dollar
amount threshold set for the Beneficiary’s Country Code as defined in Table
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBOSTOCT
Rule Name: Beneficiary with Single Transaction above Country Threshold
Customer Sub Customer Key
Priority: 7
Risk Weighting: 7
55
44. WBOSTOCT: Originator with Single Transaction above Country
Threshold
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-11-O and G4-11-O. This scenario detects if an Originator has a transaction
with a USD amount greater than the single dollar amount threshold set for the Originator’s Country.
When transactions for an Originator meet these criteria, an Alert will be raised.
For example, in the Look Back Project a transaction threshold was derived for the UAE as a single
dollar amount threshold of USD 2,100,000.00. If an Originator in the UAE had a single transaction with
a USD amount greater than USD 2,100,000.00 the Alert would generate.
This Rule is designed for those entities (i.e. Originators) that have only a single transaction in a month,
and therefore are unlikely to hit the volume thresholds embedded within other Rules. This Rule will
capture large dollar amount transactions, and subject them to a review.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Originator Country is identified by the following logic:
Transaction-daily_yyyymmdd.originator_country_code (If Blank, check next)
• The Originator Single Transaction Dollar Amount Threshold for each Country Code is
defined in Table
• If the Originator field is blank, do not run the rule.
• If the Originator field is populated but the Originator Country field is blank, default to the
Country Code = NQ representing Null.
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
Pseudo Code
IF
An Originator has a transaction in a month with a USD amount greater than the single dollar
amount threshold set for the Originator’s Country Code as defined in Table
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBOSTOCT
Rule Name: Originator with Single Transaction above Country Threshold
Customer Sub Customer Key
Priority: 7
Risk Weighting: 7
56
45. WBOHRBHR: Originator in High Risk Country/Jurisdiction with
Transactions above Country Thresholds with Beneficiaries in High
Risk
Countries/Jurisdictions
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-8-O and G4-8-O.
This scenario detects if an Originator in a High Risk Country or Jurisdiction transacts with Beneficiaries
in High Risk Countries or Jurisdictions (excluding the Originator’s Country) , and if so, whether the total
dollar amount of the Originator’s transactions is greater than the total dollar amount threshold set for
the Originator’s Country in a month. When transactions for an Originator meet these criteria, an Alert
will be raised.
For example, if an Originator in the UAE (which is a High Risk Country) transacts with Beneficiaries in
Armenia for a total amount of USD 1,000,000.00, and transacts with Beneficiaries in Bangladesh for a
total amount of USD 4,000,000.00, both of which are High Risk, the transactions in a month aggregate
to a total dollar amount of USD 5,000,000.00 which exceeds the Originator’s threshold dollar amount
for UAE (USD 3,200,000.00) the alert would generate.
Originator Beneficiaries
UAE USD 3,200,000.00 Armenia USD 1,000,000.00
(Threshold) UAE USD 2,000,000.00 excluded
Bangladesh USD 4,000,000.00
USD 5,000,000.00
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Originator Country is identified by the following logic:
Transaction-daily_yyyymmdd.originator_country_code
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Beneficiary Country is identified by the following logic:
Transaction-daily_yyyymmdd.beneficiary_country_code
• The Originator Dollar Amount Threshold for each Country Code is defined in Table
• If the Originator field is blank, do not run the rule.
• The Originator’s Home Country is defined as a Beneficiary Country Code identical to the
Originator Country Code, e.g. Beneficiary Country Code = AE and Originator Country Code =
AE
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
57
Pseudo Code
IF
The Originator Country Code has the Internal Risk Flag = Yes
AND
The Beneficiary Country Code has the Internal Risk Flag = Yes
AND
The Originator has only one country code not including blank or null and the total USD dollar
amount of the total number of transactions for an Originator in a High Risk Country with all
Beneficiaries in High Risk Countries (excluding the Originator’s Home Country) in a month is
greater than the total dollar amount threshold set for the Originator’s Country Code as defined
in the Country Threshold Table.
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBOHRBHR
Rule Name: Originator in High Risk Country with Transactions above Country Thresholds with
Beneficiary in High Risk Country
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
58
46. WBBHROHR: Beneficiary in High Risk Country/Jurisdiction with
Transactions above Country Thresholds from Originators in High Risk
Countries/Jurisdictions
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Look Back Project Rule G3-8-B and G4-8-B.
This scenario detects if an Beneficiary in a High Risk Country or Jurisdiction transacts with Originators
in High Risk Countries or Jurisdictions (excluding the Beneficiary’s Country) , and if so, whether the
total dollar amount of the Beneficiary’s transactions is greater than the total dollar amount threshold
set for the Beneficiary’s Country in a month. When transactions for a Beneficiary meet these criteria,
an Alert will be raised.
For example, if a Beneficiary in the UAE (which is a High Risk Country) receives payments from
Originators in Armenia for a total amount of USD 1,000,000.00, and receives payments from
Originators in Bangladesh for a total amount of USD 4,000,000.00, both of which are High Risk, the
transactions in a month aggregate to a total dollar amount of USD 5,000,000.00 which exceeds the
Beneficiary’s threshold dollar amount of USD 2,300,000.00 the Alert would generate.
Originator Beneficiaries
UAE USD 2,300,000.00 Armenia USD 1,000,000.00
(Threshold) UAE USD 2,000,000.00 excluded
Bangladesh USD 4,000,000.00
USD 5,000,000.00
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at Sub Customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• The Originator Country is identified by the following logic:
Transaction-daily_yyyymmdd.originator_country_code
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Beneficiary Country is identified by the following logic:
Transaction-daily_yyyymmdd.beneficiary_country_code
• The Beneficiary Dollar Amount Threshold for each Country Code is defined in Table
• If the Beneficiary field is blank, do not run the rule.
• The Beneficiary’s Home Country is defined as a Beneficiary Country Code identical to the
Originator Country Code, e.g. Beneficiary Country Code = AE and Originator Country Code =
AE
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
59
Pseudo Code
IF
The Beneficiary Country Code has the Internal Risk Flag = Yes
AND
The Originator Country Code has the Internal Risk Flag = Yes
AND
The Beneficiary’s Country Code does not equal their Home Country
AND
The Originator’s Country Code does not equal their Home Country
AND
The Originator has only one country code not including blank or null and the total USD dollar
amount of the total number of transactions for a Beneficiary in a High Risk Country from all
Originators in High Risk Countries (excluding the Beneficiary’s Home Country) in a month is
greater than the total dollar amount threshold set for the Beneficiary’s Country Code as
defined in the Country Threshold Table.
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBBHROHR
Rule Name: Beneficiary in High Risk Country with Transactions above Country Thresholds
from an Originator in High Risk Country
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
60
47. WBCMRTOB: Roundtrips between an Originator and
Beneficiary
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank and corresponds to Look
Back Project Rules G3-13 and G4-13.
The scenario is run monthly to detect when the transactions with an amount greater than USD 50,000
from an Originator to a Beneficiary are matched to the transactions from that Beneficiary back to the
Originator within the same month. Transactions are not matched on amount. When an originator
transacts with a beneficiary and then the beneficiary transacts back with the originator an alert will be
raised in the case management software. The alert will have a code of ‘WBCMRTOB’ and a name of
‘Round Trips between an Originator and Beneficiary’. The priority displayed on the end User screen
will be set to 5, to indicate that the alert is Medium priority.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at account level.
Assumptions
• Originator is defined as one of the following:
TRANSACTION_DAILY_YYYYMMDD.originator (if blank, check next)
TRANSACTION_DAILY_YYYYMMDD.ordering_bank (if blank, check next)
TRANSACTION_DAILY_YYYYMMDD.sending_bank (if blank, check next)
TRANSACTION_DAILY_YYYYMMDD.debit_party
• Beneficiary is defined as one of the following:
TRANSACTION_DAILY-YYYYMMDD.beneficiary (if blank, check next)
TRANSACTION_DAILY_YYYYMMDD.beneficiary_bank (if blank, check next)
TRANSACTION_DAILY_YYYYMMDD.credit_party
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Transaction 1 amount > USD 50,000.00
AND
Transaction 2 amount > USD 50,000.00
AND
Transaction 1 originator = Transaction 2 beneficiary
AND
Transaction 1 beneficiary = Transaction 2 originator
AND
Organisation Unit = ‘US WB’
THEN
Create Account Alert
Rule Code WBCMRTOB,
Rule Name Round Trips between an Originator and Beneficiary,
Customer Customer Key,
61
48. WBOST6MA: Originator Structuring – 6 Month Aggregation
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Assist//ck 6 Month Structuring ORG Rule.
This scenario detects if the number of transactions for an Originator with a USD dollar amount of
between USD8,975.00 and USD 10,000.00 is equal to/greater than 10 within the last 180 days or 6
months with a total aggregate amount equal to/greater than USD 100,000.00. When transactions for
an Originator meet these criteria, an Alert will be raised. This rule is designed to complement the
existing rule WBCMORST which aggregates transactions within the last 30 days, and the new 2 Month
Aggregation rule. An Originator may not have enough transactions meeting the rule criteria within any
30 or 60 day period, but may have enough transactions meeting the rule thresholds when the
aggregation is over 6 months.
Rule No. Rule Name Aggregation
Period
Dollar Amount No.
Transactions
Aggregation
Amount
USG304O ORG
Structuring
30 days 8975.00 –
10,000.00
10
ORG
Structuring – 2
month
aggregation
60 days 8975.00 –
10,000.00
10 100,000.00
ORG
Structuring – 6
month
aggregation
180 days 8975.00 –
10,000.00
10 100,000.00
The typical 6 month scenario would be the following:
Scenario
Originator (ABC Inc.) has the following number of transactions meeting the Originator Structuring
criteria in each month.
Value Date Amount Originator
04/11/2006 9,985.00 ABC Inc
04/13/2006 9,985.00 ABC Inc
04/18/2006 9,985.00 ABC Inc.
05/25/2006 9,985.00 ABC Inc
05/30/2006 9,985.00 ABC Inc
06/05/2006 9,985.00 ABC Inc
07/22/2006 9,985.00 ABC Inc
08/22/2006 9,985.00 ABC Inc.
08/22/2006 9,985.00 ABC Inc
09/22/2006 9,985.00 ABC Inc
09/29/2006 9,985.00 ABC Inc
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
62
The data will be monitored at sub-customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• If the Originator field is blank, do not run the rule.
* If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
• While the rule will be run monthly, it will aggregate transactions for the preceding 180 days or
6 previous months. For example, if the monthly aggregations are run on 2 October 2006, the
transactions to be aggregated are the following:
September 2006
August 2006
July 2006
June 2006
May 2006
April 2006
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for an Originator is equal to/greater than 10 within the last 180 day
period or 6 preceding months
AND
The USD dollar amount of each of the 10 transactions is between USD 8,975.00 and USD
10,000.00
AND
The total USD dollar amount of the total number of such transactions for a Beneficiary is equal
to/greater than USD 100,000.00
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBOST6MA
Rule Name: Originator Structuring – 6 Month Aggregation
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
63
49. WBBST6MA: Beneficiary Structuring – 6 Month Aggregation
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Assist//ck 6 Month Structuring 3rd Party Rule.
This scenario detects if the number of transactions for a Beneficiary with a USD dollar amount of
between USD 8,975.00 and USD 10,000.00 is equal to/greater than 10 within the last 180 days or 6
months with a total aggregate amount equal to/greater than USD 100,000.00. When transactions for a
Beneficiary meet these criteria, an Alert will be raised. This rule is designed to complement the existing
rule WBCMBEST which aggregates transactions within the last 30 days, and the new 2 Month
Aggregation rule. A Beneficiary may not have enough transactions meeting the rule criteria within any
30 or 60 day period, but may have enough transactions meeting the rule thresholds when the
aggregation is over 6 months.
Rule No. Rule Name Aggregation
Period
Dollar Amount No.
Transactions
Aggregation
Amount
USG304B BEN
Structuring
30 days 8975.00 –
10,000.00
10
BEN
Structuring – 2
month
aggregation
60 days 8975.00 –
10,000.00
10 100,000.00
BEN
Structuring – 6
month
aggregation
180 days 8975.00 –
10,000.00
10 100,000.00
The typical 6 month scenario would be the following:
Scenario
Originator (ABC Inc.) has the following number of transactions meeting the Originator Structuring
criteria in each month.
Value Date Amount Beneficiary
04/11/2006 9,985.00 ABC Inc
04/13/2006 9,985.00 ABC Inc
04/18/2006 9,985.00 ABC Inc.
05/25/2006 9,985.00 ABC Inc
05/30/2006 9,985.00 ABC Inc
06/05/2006 9,985.00 ABC Inc
07/22/2006 9,985.00 ABC Inc
08/22/2006 9,985.00 ABC Inc.
08/22/2006 9,985.00 ABC Inc
09/22/2006 9,985.00 ABC Inc
09/29/2006 9,985.00 ABC Inc
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at sub-customer level.
64
Assumptions
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• If the Beneficiary field is blank, do not run the rule.
* If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
• While the rule will be run monthly, it will aggregate transactions for the preceding 180 days or
6 previous months. For example, if the monthly aggregations are run on 2 October 2006, the
transactions to be aggregated are the following:
September 2006
August 2006
July 2006
June 2006
May 2006
April 2006
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for a Beneficiary is equal to/greater than 10 within the last 180 day
period or 6 preceding months
AND
The USD dollar amount of each of the 10 transactions is between USD 8,975.00 and USD
10,000.00
AND
The total USD dollar amount of the total number of such transactions for a Beneficiary is equal
to/greater than USD 100,000.00
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBOST6MA
Rule Name: Beneficiary Structuring – 6 Month Aggregation
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
65
50. WBSOSB6MA: Multiple Transactions - Same Originator/Same
Beneficiary – 6 Month Aggregation
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Assist//ck 6 Month Beneficiary: Originator Risk 3 Rule.
This scenario detects multiple transactions with the same Originator and the same Beneficiary. For
example:
Originator Dell Computer Corporation
Beneficiary Dell Computer Corporation
If the number of transactions for the same Originator and the same Beneficiary with a USD dollar
amount equal to/greater than USD 2,999.99 is equal to/greater than 10 within the last 180 days or 6
months with a total aggregate amount equal to/greater than USD 100,000.00, an Alert will be raised.
This rule is designed to complement the existing rule WBCMSOSB which aggregates transactions
within the last 30 days, and the new 2 Month Aggregation rule. An Originator may not have enough
transactions meeting the rule criteria within any 30 or 60 day period, but may have enough
transactions meeting the rule criteria when the aggregation is over 6 months.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at sub-customer level.
Assumptions
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
• While the rule will be run monthly, it will aggregate transactions for the preceding 180 days or
6 previous months.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Number of same Originator and same Beneficiary transactions is equal to/greater than 10
within the last 180 day period or 6 preceding months
AND
The USD dollar amount of each of the 10 transactions is equal to/greater than USD 2,999.99
AND
The total USD dollar amount of the total number of such transactions between an Originator
and a Beneficiary is equal to/greater than USD 100,000.00
AND
Bank to Bank Flag = N or Blank
AND
Organisation Unit = US WB
66
THEN
Create Sub Customer Alert
Rule Code WBSOSB6MA
Rule Name Multiple Transactions – Same Originator/Same Beneficiary – 6 Month Aggregation
Customer Sub Customer Key
Priority 10
Risk Weighting 7
67
51. WBOMT2MA: Originator with Multiple Transactions – 2 Month
Aggregation
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Assist//ck 2 Month Originator Risk 3 Rule. This scenario detects if the number of transactions for an
Originator with a USD dollar amount equal to/greater than USD 2,999.99 is equal to/greater than 5
within the last 60 days or 2 months with a total aggregate amount equal to/greater than USD
50,000.00. When transactions for an Originator meet these criteria, an Alert will be raised. This rule is
designed to complement the existing rules which aggregate transactions within the last 30 days. An
Originator may not have enough transactions meeting those rule criteria within any 30 day period, but
may have enough transactions meeting the rule thresholds when the aggregation is over 2 months.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at sub-customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• If the Originator field is blank, do not run the rule
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
• While the rule will be run monthly, it will aggregate transactions for the preceding 60 days or 2
previous months.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for an Originator is equal to/greater than 5 within the last 60 day
period or 2 preceding months
AND
The USD dollar amount of each of the 5 transactions is equal to/greater than USD 2,999.99
AND
The total USD dollar amount of the total number of such transactions for an Originator is equal
to/greater than USD 50,000.00
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBOMT2MA
Rule Name: Originator with Multiple Transactions – 2 Month Aggregation
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
68
52. WBOMT6MA: Originator with Multiple Transactions – 6 Month
Aggregation
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Assist//ck 6-Month Originator Risk 3 Rule. This scenario detects if the number of transactions for an
Originator with a USD dollar amount equal to/greater than USD 2,999.99 is equal to/greater than 10
within the last 180 days or 6 months with a total aggregate amount equal to/greater than USD
100,000.00. When transactions for an Originator meet these criteria, an Alert will be raised.
This rule is designed to complement the existing rules which aggregate transactions within the last 30
days. An Originator may not have enough transactions meeting those rule criteria within any 30 day
period, but may have enough transactions meeting the rule thresholds when the aggregation is over 6
months.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at sub-customer level.
Assumptions
• The Originator is identified by the following logic: Transaction_daily_yyyymmdd.originator
• If the Originator field is blank, do not run the rule
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
• While the rule will be run monthly, it will aggregate transactions for the preceding 180 days or
6 previous months.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for an Originator is equal to/greater than 10 within the last 180 day
period or 6 preceding months
AND
The USD dollar amount of each of the 10 transactions is equal to/greater than USD 2,999.99
AND
The total USD dollar amount of the total number of such transactions for an Originator is equal
to/greater than USD 100,000.00
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBOMT6MA
Rule Name: Originator with Multiple Transactions – 6 Month Aggregation
69
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
70
53. WBBMT2MA: Beneficiary with Multiple Transactions – 2 Month
Aggregation
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Assist//ck 2-month Beneficiary Risk 3 Rule.
This scenario detects if the number of transactions for a Beneficiary with a USD dollar amount equal
to/greater than USD 2,999.99 is equal to/greater than 5 within the last 60 days or 2 months with a total
aggregate amount equal to/greater than USD 50,000.00. When transactions for a Beneficiary meet
these criteria, an Alert will be raised.
This rule is designed to complement the existing rules which aggregate transactions within the last 30
days. A Beneficiary may not have enough transactions meeting those rule criteria within any 30 day
period, but may have enough transactions meeting the rule thresholds when the aggregation is over 2
months.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at sub-customer level.
Assumptions
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• If the Beneficiary field is blank, do not run the rule
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
• While the rule will be run monthly, it will aggregate transactions for the preceding 60 days or 2
previous months.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for a Beneficiary is equal to/greater than 5 within the last 60 day
period or 2 preceding months
AND
The USD dollar amount of each of the 5 transactions is equal to/greater than USD 2,999.99
AND
The total USD dollar amount of the total number of such transactions for an Originator is equal
to/greater than USD 50,000.00
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBBMT2MA
Rule Name: Beneficiary with Multiple Transactions – 2 Month Aggregation
71
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
72
54. WBBMT6MA: Beneficiary with Multiple Transactions – 6 Month
Aggregation
High Level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank, and corresponds to the
Assist//ck 6-Month Beneficiary Risk 3 Rule. This scenario detects if the number of transactions for a
Beneficiary with a USD dollar amount equal to/greater than USD 2,999.99 is equal to/greater than 10
within the last 180 days or 6 months with a total aggregate amount equal to/greater than USD
100,000.00. When transactions for a Beneficiary meet these criteria, an Alert will be raised.
This rule is designed to complement the existing rules which aggregate transactions within the last 30
days. A Beneficiary may not have enough transactions meeting those rule criteria within any 30 day
period, but may have enough transactions meeting the rule thresholds when the aggregation is over 6
months.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at sub-customer level.
Assumptions
• The Beneficiary is identified by the following logic: Transaction_daily_yyyymmdd.beneficiary
• If the Beneficiary field is blank, do not run the rule
• If the TXN_Source_Type_Code on the transaction equals either 014, or 016, or 514
representing a Book Transfer then count only a single occurrence of the Transaction
Reference Number (TRN No.) for either the transaction threshold or the dollar amount
threshold.
• While the rule will be run monthly, it will aggregate transactions for the preceding 180 days or
6 previous months.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
The number of transactions for a Beneficiary is equal to/greater than 10 within the last 180
days period or 6 preceding months
AND
The USD dollar amount of each of the 10 transactions is equal to/greater than USD 2,999.99
AND
The total USD dollar amount of the total number of such transactions for a Beneficiary is equal
to/greater than USD 100,000.00
AND
Bank to Bank Flag = N or Blank
AND
Organization Unit = US WB
THEN
Create Sub Customer Alert
Rule Code: WBBMT6MA
Rule Name: Originator with Multiple Transactions – 6 Month Aggregation
Customer Sub Customer Key
Priority: 8
Risk Weighting: 7
73
55. SGTRVRCU: Transactions Beyond x% Variance
High Level Description
This scenario detects the variance between the actual transactions of a customer and the customer‘s
transactional activity based on past twelve months.
When the actual monthly transaction amount for a customer across all his accounts is greater than or
equal to a defined value and it is greater than a defined average monthly amount based on past 12
months for the customer who has a relationship for 12 months, an Alert will be raised.
Execution Frequency
This detection scenario will be run on a monthly basis.
Aggregation Level
The data will be monitored at customer level.
Assumptions
Anticipated Volume is loaded into system for generation of alerts based on expected volume. If not,
only the second condition would be executed.
Pseudo Code
IF
(Customer Acquisition date > 365 days
AND
Actual Monthly Txn Amt >= Threshold value
AND Actual Monthly Txn Amt is >= Factor * Avg monthly amt based on past 12 months
AND
Organisation Unit = ”WB/CB/SA/PB‘
THEN
Create Customer Alert
Rule Code SGTRVRCU,
Rule Name Transactions Beyond x% Variance at Customer Level
74
56. USG01THC: Transactions above Threshold – Monthly
High Level Description
This scenario is run for the Standard Chartered Bank’s GB Wholesale Bank and detects the variance
between the actual transactions in a customer’s account for a designated interval (e.g. monthly) and
the customer’s expected level of transactional activity as reported by the Relationship Manager
through the KYC due diligence process on the eKYC Form and US Addendum.
The Rule will compare the customer’s monthly Volume and Amount to the eKYC Details (Transactional
Profile). If there is a variance or breach of the eKYC Details, the Rule will then match the eKYC Details
Volume or Amount to a Tier/Band in the Volume Group Table to determine the Percent tolerance by
which the variance can exceed the eKYC Details. If the variance is greater than the result of
calculating the Percent tolerance and adding it to the eKYC value, an Alert will be raised.
Execution Frequency
This detection scenario will be run on a monthly basis.
Aggregation Level
The data will be monitored at account level.
Assumptions
• eKYC Details
New columns will be added to the Customer and Account records to show expected monthly
credit volume, expected monthly debit volume, expected monthly USD credits and expected
monthly USD debits.
• Volume Groups Table
A Volume Group Table is to be added to the Norkom Alchemist database.
The logical relationship between each of the eKYC Details fields and the Volume Group Table
Types is
As follows:
eKYC Details Volume Group Table
Expected Monthly Credit Volume Column A Type = TXN
Expected Monthly Debit Volume Column A Type = TXN
Expected Monthly Base Currency Credits Column A Type = Amount
Expected Monthly Base Currency Debits Column A Type = Amount
• New Account
A new account is defined as a DDA account which has been open less than or equal to 180
days. (Current System Date less Account Open Date <= 180 days)
• Regular Account
A regular account is defined as any DDA account which has been opened more than 180
days. (Current System Date less Account Open Date >= 180 days)
Pseudo Code
IF
The date account opened is greater than 180 days from current date (defined as not a new
account)
AND
Total Number of Current Month Credits = > the result of calculating the Percent Tolerance
(Column E) set for the Account’s eKYC Profile Expected Monthly Credits Volume (Column A
75
Type TXN) Tier/Band (Column C to D) in the Volume Group Table plus the eKYC value. [For
example, .75 x 100 = 75 + 100 = 175]
OR
Total Number of Current Month Debits = > the result of calculating the Percent Tolerance
(Column E) set for the Account’s eKYC Profile Expected Monthly Debits Volume (Column A
Type TXN) Tier/Band (Column C to D) in the Volume Group Table plus the eKYC value. [For
example, .75 x 100 = 75 + 100 = 175]
OR
Total Amount of Current Month USD Credits = > the result of calculating the Percent
Tolerance (Column E) set for the Account’s eKYC Profile Expected Monthly USD Credits
(Column A Type Amount). Tier/Band (Column C to D) in the Volume Group Table plus the
eKYC value. [For example .50 x 1,000,000.00 = 500,000.00 + 1,000,000.00 = 1,500,000.00]
OR
Total Amount of Current Month USD Debits = > the result of calculating the Percent Tolerance
(Column E) set for the Account’s eKYC Profile Expected Monthly USD Debits (Column A Type
Amount) Tier/Band (Column C to D) in the Volume Group Table plus the eKYC value. [For
example .50 x 1,000,000.00 = 500,000.00 + 1,000,000.00 = 1,500,000.00]
AND
Organisation Unit = ‘US WB’
THEN
Create Account Alert
Rule Cod: USG01THC,
Rule Name: Transactions Above Threshold - Monthly
Customer Key,
Account Key,
Priority: 5
Risk Weighting: 6
76
57. PBNEWCUS: New Customer Monitoring
High Level Definition
This scenario monitors New to Bank customers within the last 6 months and an alert is generated
when monthly or debit turnover is greater than or equal to Threshold amount.
Execution Frequency
This detection scenario will be executed on a Monthly basis.
Aggregation Level
The data will be monitored at the Customer Level
Assumptions
New to Bank Customer in the last 6 months
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Customer Acquisition date < 6 months
AND
Total Debit in a Month >= Threshold
OR
Total Credit in a Month >= Threshold
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Customer Alert
Rule Code: PBNEWCUS
Rule Name: New Customer Monitoring
77
58. GRWASHTR: Wash Transactions at a Customer Level
High Level Description
This scenario monitors accounts where a transaction credited to an account has been removed with a
variance of 10% in the past two months
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at customer level.
Assumptions
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Customer Acquisition date > 60 days
AND
Total Debit over 2 months >= Threshold Value
AND
Total Credit over 2 months is BETWEEN ±10% * Debit Amount
AND
Actual Monthly Txn Amt is >= Factor * Avg monthly amt based on past 2 months
AND
Organization Unit =”WB/CB/SA/PB”
THEN
Create Account Alert
Rule Code GRWASHTR,
Rule Name Transactions Beyond 10% Variance at Customer Level
78
59. CBCASTXN: Cash Transaction Monitoring
High Level Description
This scenario detects customers who have cash deposits or withdrawals for at least a day greater than
or equal to the Threshold Amount
Execution Frequency
This detection scenario will be run on a daily basis
Aggregation Level
The data will be monitored at customer level.
Assumptions
• Cash transaction types need to be populated in aggregation parameters table.
• Thresholds need to be updated in the Scenario Values table.
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
There is at least one instance in a day
AND
Total Debit over 1 day >= Threshold Value
OR
Total Credit over 1 day >= Threshold Value
AND
Organization Unit =”WB/CB/SA/PB”
THEN
Create Account Alert
Rule Code CBCASTXN,
Rule Name Cash Transaction Monitoring
79
60. GLCASBLR: Boiler Room Scam for Cash
High Level Description
This scenario detects accounts whose country of incorporation is classified as high risk, the
incorporation date is less than 12 months, the count and amount of monthly incoming TT transactions
Is >= defined and Threshold1 and amount of monthly debit transactions is >= Threshold2 and the
business is Investment or Trading. An alert is generated in this case.
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at customer level.
Assumptions
• List of Country of Incorporation that is High Risk would be provided.
• ITT transactions are also mapped to TT0019 AML Code
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Country of Incorporation is High Risk AND ISIC Code is within defined
AND
Account since Opened Date < 12 months
AND
Count of monthly Credit TT
Transactions is >= defined
AND
Amount of monthly Credit TT transactions >= Threshold1
AND
Amount of monthly Debit Transactions >= Threshold2
AND
Organization Unit =”WB/CB/SA/PB‘
THEN
Create Account Alert
Rule Code: GLCASBLR
Rule Name: Boiler Room Scam for Cash
80
61. INIBACSA: Cash Deposit Structuring A
High Level Description
This scenario detects whether all the accounts of a customer witness a cash deposit between
threshold 1 and threshold 2 and the aggregated amount is greater than aggregated threshold
This scenario was made specifically for cases where cash deposits just below INR 10,00,000 were
split into multiple accounts in a month & frequent cash deposits just below INR 10,00,000.
Alert would be raised when
1. Cash deposits in amounts ranging between INR 9,00,000/- to INR 9,99,999.99) in multiple accounts
of the customer greater than [N] times in a month
2. Cash transactions (deposits) in amounts ranging between INR 9,00,000/- to INR 9,99,999.99)
greater than [N] times in [Y] days
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at customer level.
Assumptions
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Amount of each cash deposit should be between threshold 1 and threshold 2 (additional
variables to be created in Scenario values table)
AND
Aggregated amount of such transactions > Aggregated threshold
AND
Number of such transactions > defined count
THEN
Create Account Alert
Rule Code: INIBACSA
Rule Name: Cash Structuring Deposit A
81
62. INIBACSB: Cash Deposit Structuring B
High Level Description
This scenario detects whether all the accounts of a customer witness a cash deposit between
threshold 1 and threshold 2 and the aggregated amount is greater than aggregated threshold
This scenario was made specifically for cases where cash deposits just below INR 50,000 were split
into multiple accounts.
Alert would be raised when
1. Deposit of cash in the account in amounts ranging between INR 40,000/- to INR 49,999/- greater
than [N] times in [Y] days
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at customer level.
Assumptions
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Amount of each cash deposit should be between threshold 1 and threshold 2 (additional
variables to be created in Scenario values table)
AND
Aggregated amount of such transactions > Aggregated threshold
AND
Number of such transactions > defined count
THEN
Create Account Alert
Rule Code: INIBACSB
Rule Name: Cash Structuring Deposit A
82
63. GRTRNEKY: New Account Monitoring for First 6 Months
High Level Description
This scenario detects the variance between the actual transactions in a new customer’s account
(during the first six months after opening) for a designated interval (e.g. monthly) and the customer’s
expected level of transactional activity as reported by the Relationship Manager through the KYC due
diligence process on the eKYC Form and US Addendum. While the Rule logic is identical to Rule
GRTRTAKYC, this Rule is intended to highlight a breach in a new account, which is subject to a
heightened level of monitoring during the first six months after account opening, as the pattern of
customer activity is identified and confirmed to the expected level as reported in the eKYC. During the
life cycle of an account, it would be monitored under Rule GRTRNEKYC for its first six months, and
then under Rule GRTRTAKYC for the remainder of the time it is open.
Original rule definition contained a set parameter variance of +/- 10%. However, Compliance has
introduced a Volume Group Table in Assist//CK to deal with the fact that a 10% variance on a
customer expected to do 10 transactions a month is very different from a 10% variance on another
customer expected to do 1,000 transactions a month. As part of the Alert management strategy, this
rule is therefore being amended to replace the fixed 10% parameter with a variable Percent parameter
from the Volume Group Table.
The Rule will now compare the customer’s monthly Volume and Amount to the eKYC Details
(Transactional Profile). If there is a variance or breach of the eKYC Details, the Rule will then match
the eKYC Details Volume or Amount to a Tier/Band in the Volume Group Table to determine the
Percent tolerance by which the variance can exceed the eKYC Details. If the variance is greater than
the result of calculating the Percent tolerance and adding it to the eKYC value, an Alert will be raised.
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at account level.
Assumptions
• eKYC Details
New columns will be added to the Customer and Account records to show expected monthly credit
volume, expected monthly debit volume, expected monthly USD credits and expected monthly USD
debits.
• Volume Groups Table
A Volume Group Table is to be added to the Norkom Alchemist database in the format attached
(Version 0.3), and to be invoked by this Rule.
The logical relationship between each of the eKYC Details fields and the Volume Group Table Types
is as follows:
eKYC Details Volume Group Table
Expected Monthly Credit Volume Column A Type = TXN
Expected Monthly Debit Volume Column A Type = TXN
Expected Monthly Base Currency Credits Column A Type = Amount
Expected Monthly Base Currency Debits Column A Type = Amount
Pseudo Code
IF
//Inclusion/Exclusion logic
Transaction Country is equal to Country Code
AND Transaction Organization Unit is equal to Organization Unit
AND Segment Code is in Segment Codes for this rule
83
AND Transaction Type is Transaction Type for this rule
AND Transaction Product is Transaction Product
//Filters
AND Pseudo Account Flag is equal to ‘N’ //Front-end: XML
AND First Leg Flag is equal to ‘Y’ //Redundancy Exclusion
AND Own Transfer Transactions Flag is equal to ‘N’ //Front-end: XML
//Monitoring period
AND Transaction Date within the preceding 6 Months
//Rule logic
AND Date account opened is less than //new account
Account Date threshold days from scenario’s run date //Front-end: XML
AND (
Monthly Credit Transaction Volume is greater than or equal to
Monthly Credit Transaction Volume Tolerance *
Expected Monthly Credit Volume //according to eKYC profile
OR Monthly Debit Transaction Volume is greater than or equal to
Monthly Debit Transaction Volume Tolerance *
Expected Monthly Debit Volume //according to eKYC profile
OR Monthly Credit Transaction Amount is greater than or equal to
Monthly Credit Transaction Amount Tolerance *
Expected Monthly Credit Amount //according to eKYC profile
OR Monthly Debit Transaction Amount is greater than or equal to
Monthly Debit Transaction Amount Tolerance *
Expected Monthly Debit Amount //according to eKYC profile
)
THEN
Create Account Alert using the following arguments:
Rule Code: GRTRTAKYC
Rule Name: Transactions above Threshold - Monthly
Customer Key
Account Key
Priority: Rule Risk Weighting * Customer Risk Rating
84
64. WBCSPDMB: Cheques Same Payee Multiple Banks
High level Description
This scenario detects if the number of checks issued by multiple DDA Accounts with the same Payee
is greater than or equal to 2 per month with an aggregate total amount of USD 25,000.00. When
checks for a Payee meet these criteria, an Alert will be raised.
For this Detection Scenario, a Check is defined as representing the following US eBBS and Norkom
Transaction Type Code.
eBBS
Code
eBBS
Description
Norkom
Code
Norkom
TXN_TYPE_DESC
010 Check Paid TT0003 Cheque
Filtering for only this Transaction Code eliminates other US eBBS transactions like electronic transfers,
service fees, system-generated interest etc. that should not count against the transaction threshold.
For this Detection Scenario, a DDA Account is defined as representing the following US eBBS Product
Code.
Code Description
100 Demand Deposit Account (001)
Filtering for only this Product Code will eliminate capturing other type accounts like Trade accounts
(400), Customer Segregated Funds accounts (500) etc.
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at sub-customer level.
Assumptions
• Alerts fired from this Detection Scenario should be discreet Alerts just for the Check rules, and
should not roll-up with rule breaches for other Detection Scenarios which may fire on the same
Payee name but which aggregate at the Transaction or KYCC level. For example:
Alert ID Alert Name Alert Details Notes
1234 Fashion House Inc. Rule: Checks - Same Payee Drawn on Multiple Banks
Rule: Checks - Same Payee Sequential Check Numbers
1235 Fashion House Inc. Rule: WBOMRDAM – Originator Multiple Round Dollar
• The Payee will be mapped in the Beneficiary field.
The Beneficiary is identified by the following logic: Transaction_Daily_yyyymmdd.Beneficiary
The DDA Account is identified by the following logic:
Transaction_Daily_yyyymmdd.Account_Source_Unique_ID
Pseudo Code
IF
//Inclusion/Exclusion logic
Transaction Country is equal to Country Code
AND Transaction Organization Unit is Organization Unit
AND Transaction Product is Transaction Product
85
AND Transaction Code is in Transaction Codes for this rule
AND Transaction Type is Transaction Type for this rule
AND Segment Code is in Segment Codes for this rule
//Filters
AND Pseudo Account Flag is equal to ‘Y’ //Front-end: XML
AND Bank to Bank Flag is equal to ‘N’ OR NULL
AND First Leg Flag is equal to ‘Y’ //( Redundancy Exclusion)
AND Beneficiary is not NULL
AND Own Transfer Transactions Flag is equal to ‘N’
//Monitoring period
AND Transaction Date within the preceding Month
//Rule Logic
AND Number of checks for same Beneficiary is greater than 1
AND Number of different Originators (DDA Accounts) is greater than 1
AND Total (Aggregate) Transaction Amount is greater than or equal to
Total (Aggregate) Transaction Amount //Front-end: XML
THEN
Create Sub Customer Alert using the following arguments:
Rule Code WBCSPDMB
Rule Name: Checks – Same Payee Multiple Banks
Sub Customer Key
Priority: Rule Risk Weighting * Customer Risk Rating
Risk Weighting: 10
86
65. WBSQCKNO – Same payee Multiple sequential cheques
High level Description
This scenario is run for the Standard Chartered Bank’s US Wholesale Bank.
This scenario filters for only Checks Paid through Customer DDA Accounts.
This scenario detects if more than 2 multiple sequentially numbered USD checks with an aggregate
total amount of USD 25,000.00 are issued from a DDA account to the same Payee within the same
month. When checks for a Payee meet these criteria, an Alert will be raised.
For this Detection Scenario, a Check is defined as representing the following US eBBS and Norkom
Transaction Type Code.
eBBS
Code
eBBS
Description
Norkom
Code
Norkom
TXN_TYPE_DESC
010 Check Paid TT0003 Cheque
Filtering for only this Transaction Code eliminates other US eBBS transactions like CHIPS, Fed Funds,
service fees, system-generated interest etc. that should not count against the transaction type.
For this Detection Scenario, a DDA Account is defined as representing the following US eBBS Product
Code.
Code Description
100 Demand Deposit Account (001)
Filtering for only this Product Code will eliminate capturing other type accounts like Trade accounts
(400), Customer Segregated Funds accounts (500) etc.
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at account level.
Assumptions
• Alerts fired from this Detection Scenario should be discreet Alerts just for the Check rules, and
should not roll-up with rule breaches for other Detection Scenarios which may fire on the same
Payee name but which aggregate at the Transaction or KYCC level. For example:
Alert ID Alert Name Alert Details Notes
1234 Fashion House Inc. Rule: Checks - Same Payee Drawn on Multiple Banks
Rule: Checks - Same Payee Sequential Check Numbers
1235 Fashion House Inc. Rule: WBOMRDAM – Originator Multiple Round Dollar
• The Payee will be mapped in the Beneficiary field. The Beneficiary is identified by the following
logic: Transaction_Daily_yyyymmdd.Beneficiary
• The DDA Account is identified by the following logic:
Transaction_Daily_yyyymmdd.Account_Source_Unique_ID
The Check Number is identified by the following logic:
Transaction_Daily_yyyywwdd.Source_Txn_Num
87
Pseudo Code
IF
//Inclusion/Exclusion logic
Transaction Country is equal to Country Code
AND Transaction Organization Unit is Organization Unit
AND Transaction Product is Transaction Product
AND Transaction Code is in Transaction Codes for this rule
AND Transaction Type is Transaction Type for this rule
AND Segment Code is in Segment Codes for this rule
//Filters
AND Pseudo Account Flag is equal to ‘N’ //Front-end: XML
AND Bank to Bank Flag is equal to ‘N’ OR NULL
AND First Leg Flag is equal to ‘Y’ //( Redundancy Exclusion)
AND Beneficiary is not NULL
AND Own Transfer Transactions Flag is equal to ‘N’
//Monitoring period
AND Transaction Date within the preceding Month
//Rule Logic
AND Number of checks for same Beneficiary is greater than 1
AND Check Numbers for that Beneficiary are sequential
AND Total (Aggregate) Transaction Amount is equal to/greater
than Total (Aggregate) Transaction Amount //Front-end: XML
THEN
Create Customer Alert using the following arguments:
Rule Code WBSQCKNO
Rule Name: Checks – Same Payee Multiple Sequential Checks
Customer Key
Account Key
Priority: Rule Risk Weighting * Customer Risk Rating
Risk Weighting: 10
88
66. GRHVLTXN: Surge in Volume Transactions
High level Description
This is a Detection Scenario which monitors surge in volume of transactions. This rule is specific to
India.
Execution Frequency
This detection scenario will be run on a monthly basis
Aggregation Level
The data will be monitored at account level.
Assumptions
The Account is identified by the following logic:
If the CUSTOMER_TYPE_CODE is __ include the transaction.
Pseudo Code
IF
//Inclusion/Exclusion logic
Transaction Country is equal to Country Code
AND Transaction Organization Unit is Organization Unit
AND Transaction Product is Transaction Product
AND Transaction Code is in Transaction Codes for this rule
AND Transaction Type is Transaction Type for this rule
AND Segment Code is in Segment Codes for this rule
//Filters
AND Pseudo Account Flag is equal to ‘Y’ //Front-end: XML
AND Bank to Bank Flag is equal to ‘N’ OR NULL
AND First Leg Flag is equal to ‘Y’ //( Redundancy Exclusion)
AND Beneficiary is not NULL
AND Own Transfer Transactions Flag is equal to ‘N’
//Monitoring period
AND Transaction Date within the preceding Month
//Rule Logic
AND Number of checks for same Beneficiary is greater than 1
AND Number of different Originators (DDA Accounts) is greater than 1
AND Total (Aggregate) Transaction Amount is greater than or equal to
Total (Aggregate) Transaction Amount //Front-end: XML
THEN
Create Sub Customer Alert using the following arguments:
Rule Code WBCSPDMB
Rule Name: Checks – Same Payee Multiple Banks
Sub Customer Key
Priority: Rule Risk Weighting * Customer Risk Rating
Risk Weighting: 10
89
67. UKLREACC: UK Lately Reactivated Accounts
High level Description
The DS rule “Lately Re-activated Accounts” has been agreed to fire on a monthly basis by the Group.
However, as per the requirement of UK Compliance, the DS rule should fire on a daily basis. This
scenario is specific to UK-WB
Execution Frequency
This detection scenario will be run on a daily basis
Aggregation Level
The data will be monitored at account level.
Assumptions
Pseudo Code
IF
Dormant Status = 2
AND
Account Reactivation Date < 1 day excluding interest
AND
Organization Unit = “GB-WB”
THEN
Create Account Alert
Rule Code: GRLREACC
Rule Name: Lately Reactivated Accounts
Account Account Key
Customer Customer Key
Priority: 8
Risk Weighting: 8
90
68. GRORGBLR: Boiler Room Scam - Originator
High Level Definition
This rule is a fine tuned version of Boiler Room Scam to enable the monitoring of transactions
irrespective of the ISIC codes. This scenario is specific to UAE.
Execution Frequency
This detection scenario will be executed on a monthly basis.
Aggregation Level
The data will be monitored at the Account Level
Assumptions
List of Country of Incorporation that is High Risk would be provided.
ITT transactions are also mapped to TT0019 AML Code
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
Number of Wire Transfer >= 3
AND
Amount of each transaction is between Threshold 1 and Threshold 2
AND
Originator Country = “US/UK/HK”
AND
Organization Unit = ”WB/CB/SA/PB‘
THEN
Create Account Alert
Rule Code: GRORGBLR
Rule Name: Boiler Room Scam - Originator
91
69. HKTALAMT – Transactions above Limit (Amount)
High Level Definition
This scenario detects if actual daily transaction amount is greater than the threshold values set for
each level of accounts.
Execution Frequency
This detection scenario will be executed on a daily basis.
Aggregation Level
The data will be monitored at the Account Level
Assumptions
The threshold values will be changed by SCB as required and as appropriate based on second level of
criteria from time to time. The second level criteria are:
1. ARM Code
2. Customer ISIC code
3. Customer ID in Hogan
4. Account Currency in Hogan
Pseudo Code
The High Level Pseudo code for deployment of this scenario is specified below:
IF
( Account_Level = L1 and txn amount aggregate > HK $ 1M
OR Account_Level = L2 and txn amount aggregate > HK $ 0.75M
OR Account_Level = L3 and txn amount aggregate > HK $ 0.5M
AND Organization Unit = “HK-WB”
AND ARM_Code NOT of WB individuals
)
OR
( Account_Level = L1 and txn amount aggregate > HK $ 0.5M
OR Account_Level = L2 and txn amount aggregate > HK $ 0.5M
OR Account_Level = L3 and txn amount aggregate > HK $ 0.25M
AND Organization Unit = “HK-WB”
AND ARM_Code of WB individuals
)
THEN
Create Account Alert
Rule Code: HKTALAMT
Rule Name: Transactions Above Limit (Amount)
Account Account Key
Customer Customer Key
Priority: 10
Risk Weighting: 8
92

  1 #!/usr/local/bin/python
  2 import psycopg2
  3 import os
  4
  5 class Upserter:
  6
  7     def __init__(self, _schemaName, _tableName, _errorFileFullPath, locally=False):
  8         self.schemaName = _schemaName
  9         self.tableName = _tableName
 10         self.errorFileFullPath = _errorFileFullPath
 11
 12         if not locally:
 13             self.conn = psycopg2.connect(database=os.environ["DATABASE"],
 14                 user=os.environ["POSTGRES_USER"],
 15                 host=os.environ["POSTGRES_HOST"],
 16                 port=os.environ["POSTGRES_PORT"],
 17                 password=os.environ["POSTGRES_PASS"])
 18         else:
 19             self.conn = psycopg2.connect(database=os.environ["DATABASE"],
 20                 user=os.environ["POSTGRES_USER"])
 21         self.conn.autocommit = True
 22
 23         self.cur = self.conn.cursor()
 24         self.discoverAllColumns()
 25         self.discoverPrimaryKeyColumns()
 26
 27     def discoverAllColumns(self):
 28         allColumns_cmd = ('SELECT attnum, attname'
 29                           ' FROM   pg_attribute'
 30                           ' WHERE  attrelid = %s::regclass'
 31                           ' AND    attnum > 0'
 32                           ' AND    NOT attisdropped'
 33                           ' ORDER  BY attnum')
 34         self.cur.execute(allColumns_cmd, (self.schemaName + '.' + self.tableName, ))
 35         allColumns_results = self.cur.fetchall()
 36         self.allColumnNames = [i[1] for i in allColumns_results]
 37
 38     def getColumnCount(self):
 39         return len(self.allColumnNames)
 40
 41     def discoverPrimaryKeyColumns(self):
 42         pkeys_cmd = ('SELECT distinct' # distinct added by me
 43                     ' pg_attribute.attname,'
 44                     ' pg_attribute.attnum,' # added
 45                     ' format_type(pg_attribute.atttypid, pg_attribute.atttypmod)'
 46                     ' FROM pg_index, pg_class, pg_attribute'
 47                     ' WHERE'
 48                     ' pg_class.oid = %s::regclass AND'
 49                     ' indrelid = pg_class.oid AND'
 50                     ' pg_attribute.attrelid = pg_class.oid AND'
 51                     ' pg_attribute.attnum = any(pg_index.indkey)'
 52                     ' AND indisprimary'
 53                     ' ORDER BY pg_attribute.attnum')
 54         self.cur.execute(pkeys_cmd, (self.schemaName + '.' + self.tableName, ))
 55         pkey_cmd_results =  self.cur.fetchall()
 56         self.pkey_col_names = [i[0] for i in pkey_cmd_results]
 57         self.pkey_col_numbers = [i[1] for i in pkey_cmd_results]
 58
 59     def checkIfPrimaryKeyExists(self, row):
 60         # Below the schema.table can't be parameterized with %s since they would then be surrounded by single-quotes.
 61         row_select_cmd = ('SELECT *'
 62                          ' FROM ' + self.schemaName + '.' + self.tableName + ''
 63                          ' WHERE'
 64                          ' ')
 65         where_clause = []
 66         execute_query_args = []
 67
 68         for idx, pkey_col_number in enumerate(self.pkey_col_numbers):
 69             where_clause.append(self.pkey_col_names[idx] + ' = %s')
 70             execute_query_args.append(row[pkey_col_number - 1])
 71         where_clause = ' AND '.join(where_clause)
 72         row_select_cmd += where_clause
 73
 74         self.cur.execute(row_select_cmd, execute_query_args)
 75         existing_rows = self.cur.fetchall()
 76         does_row_match_primary_key = False
 77         is_row_identical = False
 78         if len(existing_rows) > 0:
 79             does_row_match_primary_key = True
 80             existing_row = list(existing_rows[0])
 81             # Note: The rows will not be treated as identical if the existing row
 82             # has kea_imported='Y' while the new row has kea_imported='N'.
 83             # (same for the phx_imported column).
 84             # This means that reimporting a file means it will be reimported to phoenix.
 85             is_row_identical = (existing_row == row)
 86
 87         return does_row_match_primary_key, is_row_identical
 88
 89
 90     # row: A list of strings, one per column.
 91     def upsert(self, row):
 92         row_matches_pkey = None
 93         upsert_type = None
 94         try:
 95             if len(row) != self.getColumnCount():
 96                 raise Exception("Incorrect column count in row. Should be {0}, but was {1}".
 97                     format(self.getColumnCount(), len(row)))
 98             row_matches_pkey, row_identical = self.checkIfPrimaryKeyExists(row)
 99             # print "row_matches_pkey: " + str(row_matches_pkey)
100             if row_matches_pkey and not row_identical:
101                 self.update(row)
102                 upsert_type = 'update'
103             elif not row_matches_pkey:
104                 self.insert(row)
105                 upsert_type = 'insert'
106             else:
107                 upsert_type = 'update'
108         except Exception, e:
109             with open(self.errorFileFullPath, 'a') as errorfile:
110                 errorfile.write("Error in row {0}\n".format(str(row)))
111                 errorfile.write("    Exception: {0}\n".format(str(e)))
112             return False, None
113
114         return True, upsert_type
115
116     def insert(self, row):
117         all_columns_commas = ', '.join(self.allColumnNames)
118         # parameters_string just looks like '%s, %s, ...'
119         parameters_string = ', '.join(['%s' for rowEntry in row])
120         insert_cmd = ('INSERT INTO ' + self.schemaName + '.' + self.tableName + ''
121                       ' VALUES'
122                       ' (' + parameters_string + ')')
123         #print 'insert_cmd: ' + insert_cmd
124         #print row
125         self.cur.execute(insert_cmd, row)
126         # TODO: check for errors or exceptions.
127
128     def update(self, row):
129         set_statement = []
130         update_query_args = []
131         for idx, colName in enumerate(self.allColumnNames):
132             set_statement.append('"{0}" = %s'.format(colName)) # row[idx]
133             update_query_args.append(row[idx])
134         set_statement = ', '.join(set_statement)
135         where_statement = []
136         for idx, pkey_col_number in enumerate(self.pkey_col_numbers):
137             where_statement.append(self.pkey_col_names[idx] + ' = %s') #row[pkey_col_number - 1]
138             update_query_args.append(row[pkey_col_number - 1])
139         where_statement = ' AND '.join(where_statement)
140         update_cmd = ('UPDATE ' + self.schemaName + '.' + self.tableName + ''
141                       ' SET ' + set_statement + ''
142                       ' WHERE ' + where_statement)
143
144         #print 'update_cmd: ' + update_cmd
145         #print update_query_args
146         self.cur.execute(update_cmd, update_query_args)
147         # TODO check for errors or exceptions
148
149         # TODO: if this fails --> log the row to a .errors file or something
150
152     # IMPORTANT: This must be called following all upsert calls
153     # to commit the changes to the database!
154     def close(self):
155         self.cur.close()
156         self.conn.close()
157
158
159 if __name__ == '__main__':
160     # upserter = Upserter("cash_aml", "test_av", "av_table2")
161     upserter = Upserter("cash_aml", "test_av", "accinfocopy", "/opt/palantir/postgres-import/scripts-tmp/upsert-scripts/test-errors.log")
162     # This should update (same data though)
163     upserter.upsert(["23033079550", "ITL", "GII", "CPR", "1", "BH"])
164     # This should insert
165     upserter.upsert(["z1", "z", "z", "z", "z", "tooLongWillBeExceptionzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz    zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"    ])
166
167     upserter.close()
upsert.py                                                         

-- @c:\palanthir\trc.sql

drop table trace_log;

create table trace_log
(
id  		 integer primary key,
msg 		 varchar2(4000),
created_by   varchar2(50) default sys_context ('userenv', 'session_user'),
created_date timestamp default systimestamp
);


create or replace package trc 
is
procedure log(p_message varchar2);
end trc;
/

create or replace package body trc 
is
procedure log(v_message varchar2)
is
 PRAGMA AUTONOMOUS_TRANSACTION;
n_id integer;
begin
select seq_trace_log.nextval into n_id from dual;

insert into trace_log(id,msg) values(n_id,v_message);
commit;
end log;
end trc;
/

Python
mouse.py
# Purpose: Automatic mouse clicks at intereval specified
import ctypes
import time
import datetime

class mouse:
	def __init__(self,_max_minute,_sleep_seconds):
		self.max_minute = _max_minute
		self.sleep_seconds = _sleep_seconds
		self.min_value = 0
		self.click()

	def click(self):
		print "Started at-" + datetime.datetime.now().strftime('%H:%M:%S')
		while self.min_value < self.max_minute:			
			ctypes.windll.user32.SetCursorPos(100, 20)
			ctypes.windll.user32.mouse_event(2, 0, 0, 0, 0) # left down
			ctypes.windll.user32.mouse_event(3, 0, 0, 0, 0) # left up
			print "clicked at-" + datetime.datetime.now().strftime('%H:%M:%S')
			time.sleep(self.sleep_seconds)
			self.min_value += 1		
		print "Ended at-" + datetime.datetime.now().strftime('%H:%M:%S')
		
# Arguments
# 1 - How many loop to run the program
# 2 - Sleep time in seconds

if __name__ =="__main__":
	mouse = mouse(140,240) # 240 sec = 4 minute, 140*4 = 560 min = 9 hour 20 min 

Backup.py
import os
import shutil
from datetime import date

class backup:
	def __init__(self):
		self.backupday = date.today()
		self.newfolder = self.backupday.strftime("%B%d")
		self.fromfolder = "c:\\work"
		self.tofolder = "c:\\backup\\" + self.newfolder
		self.backup()
	
	def backup(self):
		shutil.copytree(self.fromfolder,self.tofolder) #copying
		shutil.make_archive(self.tofolder,'zip',self.tofolder) #zipping
		shutil.rmtree(self.tofolder) #cleaning

if __name__ == "__main__":
	backup = backup()

clock
# clock
# Author : Manoj Koyadan
# Date   : Jan-18-2016
# Purpose: Clock

#!/usr/bin/python
import Tkinter as tk 
from math import sin,cos,pi
import time

class clock(tk.Tk):
	def __init__(self, *args, **kwargs):
		tk.Tk.__init__(self, *args, **kwargs)
		self.title("Clock")
		self.w = tk.Canvas(self, width=520, height=520, bg="Cyan", relief= "sunken", border=10)
		self.xcentre = 270
		self.ycentre = 270
		self.radius  = 250
		self.w.pack()
		self.w.create_oval(20,20,520,520) #Main circle
		self.w.pack()
		self.w.create_oval(265,265,275,275,fill = "Black") #Centre Dot
		self.w.pack()
		
		#Second Hand
		self.w.create_line(0,0,0,0,fill = "Blue", width=1,tags="seconds")
		#Minute Hand
		self.w.create_line(0,0,0,0,fill = "Blue", width=2,tags="minute")
		#Hour Hand
		self.w.create_line(0,0,0,0,fill = "Blue", width=4,tags="hour")
		
		for i in range (0,12):
			# Drawing hour co-orninates
			self.degree = i*30
			self.w.create_line((self.radius * cos((self.degree*pi)/180)) + 270,(self.radius * sin((self.degree*pi)/180)) + 270, (240 * cos((self.degree*pi)/180)) + 270, (240 * sin((self.degree*pi)/180)) + 270,fill = "Red", width=6
		) # 
					
		for i in range (0,60):
			# Drawing minute co-orninates
			self.degree = i*6
			self.w.create_line((self.radius * cos((self.degree*pi)/180)) + 270,(self.radius * sin((self.degree*pi)/180)) + 270, (240 * cos((self.degree*pi)/180)) + 270, (240 * sin((self.degree*pi)/180)) + 270
		) # 
		
		self.change_clock()
	
	def change_clock(self):
		hour = time.localtime()[3]
		minute = time.localtime()[4]
		seconds = time.localtime()[5]
		
		# seconds
		sec_degree = seconds*6 - 90
		sec_angle = (sec_degree*pi)/180 
		sec_x = 270 + 230 * cos(sec_angle)
		sec_y = 270 + 230 * sin(sec_angle) 
		self.w.coords("seconds", (self.xcentre,self.xcentre,sec_x,sec_y))
		
		# minute
		min_degree = minute*6 - 90
		min_angle = (min_degree*pi)/180 
		min_x = 270 + 200 * cos(min_angle)
		min_y = 270 + 200 * sin(min_angle)
		self.w.coords("minute", (self.xcentre,self.xcentre,min_x,min_y))
		
		#hour
		hour_degree = hour*30 - 75
		hour_angle = (hour_degree*pi)/180 
		hour_x = 270 + 170 * cos(hour_angle) 
		hour_y = 270 + 170 * sin(hour_angle)
		self.w.coords("hour", (self.xcentre,self.xcentre,hour_x,hour_y))
		self.after(1000, self.change_clock)
			
clock = clock()
clock.mainloop()	

1 / Palaniappan/Rajesh Kannan/Sriram Sekar/Vinoth S/
Nagaraju kosripati/venkatesh Prasanna
2 / Mayavan/Arun Prasanth/Raja R/Prabhakaran Thangasamy/
Mahesh Kumar/Poornima/Moulish/Ayan

Yamini/Paul Marla/Nataraj/vidya Durai/Rajkumar Saktivelu/Sreenivasalu Keasara
VV Surya/Yanadi/Surya Shankar/Angeline/Manivasagam Subramaniyan/Sugumar Selvarasau
Fazil Mohammed/Manam Nagarjuna/Srividy Satish/Govandhan peruri/Khaja S/ 
Lavanya/Selvarasan Rangasamy/Selvam Deivanayagam/Prasanna Kumar/Anbuselvan S  / 
Rajesh Kumar Palani/Suresh Periasamy/Suchithra Hariprasad/Bhuvaneswari Rajendiran / 
Geetha K/Praveen P/Ganesh Kumar/Senthil Kumar Anbu/Srinivasa Gautham V/Sathish Kumar S/ 
Dhivya Jambukesan/Abdul Rahman M/Babu E/Sneha Vemula/Sudhakar S/Karthikeyan R D/ 
Bhuvaneshwari Govindaraj(PIP)/Vinoth Gautham(PIP)/Arunraj Muniraj/Priya Rajendran/ 
Balamurugan Chidambaram 

select total_bytes/(1024 * 1024 ) SPACE_IN_MB
from
(
select 
sum(bytes) as total_bytes
from user_segments
where SEGMENT_TYPE in ('TABLE','LOBINDEX','INDEX','LOBSEGMENT')
 and ( SEGMENT_NAME like 'FKTA_%' 
       or
      SEGMENT_NAME in ('SYS_LOB0066875925C00012$$','SYS_LOB0066875925C00013$$',
                        'SYS_LOB0062594697C00007$$')                 
       )
);

	