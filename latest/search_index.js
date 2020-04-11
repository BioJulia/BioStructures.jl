var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#BioStructures.jl-1",
    "page": "Home",
    "title": "BioStructures.jl",
    "category": "section",
    "text": "Latest Release:(Image: Latest Release) (Image: License) (Image: Documentation) (Image: BioJulia maintainer: jgreener64)Development status:(Image: Build Status) (Image: Build status) (Image: codecov.io) (Image: Documentation)"
},

{
    "location": "index.html#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "BioStructures provides functionality to read, write and manipulate macromolecular structures, in particular proteins. Protein Data Bank (PDB), mmCIF and MMTF format files can be read in to a hierarchical data structure. Spatial calculations and functions to access the PDB are also provided. It compares favourably in terms of performance to other PDB parsers - see some benchmarks online."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install BioStructures from the Julia package REPL, which can be accessed by pressing ]:add BioStructures"
},

{
    "location": "index.html#Contributing-and-questions-1",
    "page": "Home",
    "title": "Contributing and questions",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.Detailed guidance for contributing to all BioJulia packages is provided at the BioJulia Contribution Documentation.If you have a question about contributing or using this package, you are encouraged to use the BioJulia Gitter or the Bio category of the Julia discourse site."
},

{
    "location": "documentation.html#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "documentation.html#BioStructures-documentation-1",
    "page": "Documentation",
    "title": "BioStructures documentation",
    "category": "section",
    "text": "The BioStructures.jl package provides functionality to manipulate macromolecular structures, and in particular to read and write Protein Data Bank (PDB), mmCIF and MMTF files. It is designed to be used for standard structural analysis tasks, as well as acting as a platform on which others can build to create more specific tools.It compares favourably in terms of performance to other PDB parsers - see some benchmarks online and the benchmark suite. The PDB, mmCIF and MMTF parsers currently read in the whole PDB without explicit errors (with a handful of known exceptions). Help can be found on individual functions using ?function_name."
},

{
    "location": "documentation.html#Basics-1",
    "page": "Documentation",
    "title": "Basics",
    "category": "section",
    "text": "To download a PDB file:using BioStructures\n\n# Stored in the current working directory by default\ndownloadpdb(\"1EN2\")To parse a PDB file into a Structure-Model-Chain-Residue-Atom framework:julia> struc = read(\"/path/to/pdb/file.pdb\", PDB)\nProteinStructure 1EN2.pdb with 1 models, 1 chains (A), 85 residues, 754 atomsmmCIF files can be read into the same data structure with read(\"/path/to/cif/file.cif\", MMCIF). If you want to read an mmCIF file into a dictionary to query yourself (e.g. to access metadata fields), use MMCIFDict:julia> mmcif_dict = MMCIFDict(\"/path/to/cif/file.cif\")\nmmCIF dictionary with 716 fields\n\njulia> mmcif_dict[\"_entity_src_nat.common_name\"]\n1-element Array{String,1}:\n \"great nettle\"A MMCIFDict can be accessed in similar ways to a standard dictionary, and if necessary the underlying dictionary of MMCIFDict d can be accessed with d.dict. Note that the values of the dictionary are always an Array{String,1}, even if only one value was read in or the data is numerical.MMTF files can be read into the same data structure with read(\"/path/to/mmtf/file.mmtf\", MMTF). The keyword argument gzip, default false, determines if the file is gzipped. In a similar manner to mmCIF dictionaries, a MMTF file can be read into a dictionary with MMTFDict. The values of the dictionary are a variety of types depending on the MMTF specification.Refer to Downloading PDB files and Reading PDB files sections for more options.The elements of struc can be accessed as follows:Command Returns Return type\nstruc[1] Model 1 Model\nstruc[1][\"A\"] Model 1, chain A Chain\nstruc[1][\'A\'] Shortcut to above if the chain ID is a single character Chain\nstruc[\"A\"] The lowest model (model 1), chain A Chain\nstruc[\"A\"][\"50\"] Model 1, chain A, residue 50 AbstractResidue\nstruc[\"A\"][50] Shortcut to above if it is not a hetero residue and the insertion code is blank AbstractResidue\nstruc[\"A\"][\"H_90\"] Model 1, chain A, hetero residue 90 AbstractResidue\nstruc[\"A\"][50][\"CA\"] Model 1, chain A, residue 50, atom name CA AbstractAtom\nstruc[\"A\"][15][\"CG\"][\'A\'] For disordered atoms, access a specific location AtomDisordered atoms are stored in a DisorderedAtom container but calls fall back to the default atom, so disorder can be ignored if you are not interested in it. Disordered residues (i.e. point mutations with different residue names) are stored in a DisorderedResidue container. The idea is that disorder will only bother you if you want it to. See the Biopython discussion for more.Properties can be retrieved as follows:Function Returns Return type\nserial Serial number of an atom Int\natomname Name of an atom String\naltlocid Alternative location ID of an atom Char\naltlocids All alternative location IDs in a DisorderedAtom Array{Char,1}\nx x coordinate of an atom Float64\ny y coordinate of an atom Float64\nz z coordinate of an atom Float64\ncoords coordinates of an atom Array{Float64,1}\noccupancy Occupancy of an atom (default is 1.0) Float64\ntempfactor Temperature factor of an atom (default is 0.0) Float64\nelement Element of an atom (default is \"  \") String\ncharge Charge of an atom (default is \"  \") String\nresidue Residue an atom belongs to Residue\nishetero true if the residue or atom is a hetero residue/atom Bool\nisdisorderedatom true if the atom is disordered Bool\npdbline PDB ATOM/HETATM record for an atom String\nresname Residue name of a residue or atom String\nresnames All residue names in a DisorderedResidue Array{String,1}\nresnumber Residue number of a residue or atom Int\nsequentialresidues Determine if the second residue follows the first in sequence Bool\ninscode Insertion code of a residue or atom Char\nresid Residue ID of an atom or residue (full=true includes chain) String\natomnames Atom names of the atoms in a residue, sorted by serial Array{String,1}\natoms Dictionary of atoms in a residue Dict{String,AbstractAtom}\nisdisorderedres true if the residue has multiple residue names Bool\ndisorderedres Access a particular residue name in a DisorderedResidue Residue\nchain Chain a residue or atom belongs to Chain\nchainid Chain ID of a chain, residue or atom String\nresids Sorted residue IDs in a chain Array{String,1}\nresidues Dictionary of residues in a chain Dict{String,AbstractResidue}\nmodel Model a chain, residue or atom belongs to Model\nmodelnumber Model number of a model, chain, residue or atom Int\nchainids Sorted chain IDs in a model or structure Array{String,1}\nchains Dictionary of chains in a model or structure Dict{String,Chain}\nstructure Structure a model, chain, residue or atom belongs to ProteinStructure\nstructurename Name of the structure an element belongs to String\nmodelnumbers Sorted model numbers in a structure Array{Int,1}\nmodels Dictionary of models in a structure Dict{Int,Model}The strip keyword argument determines whether surrounding whitespace is stripped for atomname, element, charge, resname and atomnames (default true).The coordinates of an atom can be set using x!, y!, z! and coords!."
},

{
    "location": "documentation.html#Manipulating-structures-1",
    "page": "Documentation",
    "title": "Manipulating structures",
    "category": "section",
    "text": "Elements can be looped over to reveal the sub-elements in the correct order:for mod in struc\n    for ch in mod\n        for res in ch\n            for at in res\n                # Do something\n            end\n        end\n    end\nendModels are ordered numerically; chains are ordered by chain ID character ordering, except the empty chain is last; residues are ordered by residue number and insertion code with hetero residues after standard residues; atoms are ordered by atom serial. If you want the first sub-element you can use first. For example first(struc[1]) gets the first chain in model 1. Since the ordering of elements is defined you can use the sort function. For example sort(res) sorts a list of residues as described above, or sort(res, by=resname) will sort them alphabetically by residue name.collect can be used to get arrays of sub-elements. collectatoms, collectresidues, collectchains and collectmodels return arrays of a particular type from a structural element or element array. Since most operations should use a single version of an atom or residue, disordered entities are not expanded by default and only one entity is present in the array. This can be changed by setting expand_disordered to true in collectatoms or collectresidues.Selectors are functions passed as additional arguments to these functions. Only elements that return true when passed to all the selector are retained. For example:Command Action Return type\ncollect(struc[\'A\'][50]) Collect the sub-elements of an element, e.g. atoms from a residue Array{AbstractAtom,1}\ncollectresidues(struc) Collect the residues of an element Array{AbstractResidue,1}\ncollectatoms(struc) Collect the atoms of an element Array{AbstractAtom,1}\ncollectatoms(struc, calphaselector) Collect the Cα atoms of an element Array{AbstractAtom,1}\ncollectatoms(struc, calphaselector, disorderselector) Collect the disordered Cα atoms of an element Array{AbstractAtom,1}The selectors available are:Function Acts on Selects for\nstandardselector AbstractAtom or AbstractResidue Atoms/residues arising from standard (ATOM) records\nheteroselector AbstractAtom or AbstractResidue Atoms/residues arising from hetero (HETATM) records\natomnameselector AbstractAtom Atoms with atom name in a given list\ncalphaselector AbstractAtom Cα atoms\ncbetaselector AbstractAtom Cβ atoms, or Cα atoms for glycine residues\nbackboneselector AbstractAtom Atoms in the protein backbone (CA, N, C and O)\nheavyatomselector AbstractAtom Non-hydrogen atoms\nhydrogenselector AbstractAtom Hydrogen atoms\nresnameselector AbstractAtom or AbstractResidue Atoms/residues with residue name in a given list\nwaterselector AbstractAtom or AbstractResidue Atoms/residues with residue name HOH\nnotwaterselector AbstractAtom or AbstractResidue Atoms/residues with residue name not HOH\ndisorderselector AbstractAtom or AbstractResidue Atoms/residues with alternative locations\nallselector AbstractAtom or AbstractResidue All atoms/residuesTo create a new atomnameselector or resnameselector:cdeltaselector(at::AbstractAtom) = atomnameselector(at, [\"CD\"])It is easy to define your own atom, residue, chain or model selectors. The below will collect all atoms with x coordinate less than 0:xselector(at) = x(at) < 0\ncollectatoms(struc, xselector)Alternatively, you can use an anonymous function:collectatoms(struc, at -> x(at) < 0)countatoms, countresidues, countchains and countmodels can be used to count elements with the same selector API. For example:julia> countatoms(struc)\n754\n\njulia> countatoms(struc, calphaselector)\n85\n\njulia> countresidues(struc, standardselector)\n85\n\njulia> countatoms(struc, expand_disordered=true)\n819The amino acid sequence of a protein can be retrieved by passing an element to AminoAcidSequence with optional residue selectors:julia> AminoAcidSequence(struc[\'A\'], standardselector)\n85aa Amino Acid Sequence:\nRCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRCThe gaps keyword argument determines whether to add gaps to the sequence based on missing residue numbers (default true). See BioSequences.jl and BioAlignments.jl for more on how to deal with sequences. For example, to see the alignment of CDK1 and CDK2 (this example also makes use of Julia\'s broadcasting):julia> struc1, struc2 = retrievepdb.([\"4YC6\", \"1HCL\"])\n2-element Array{ProteinStructure,1}:\n ProteinStructure 4YC6.pdb with 1 models, 8 chains (A,B,C,D,E,F,G,H), 1420 residues, 12271 atoms\n ProteinStructure 1HCL.pdb with 1 models, 1 chains (A), 294 residues, 2546 atoms\n\njulia> seq1, seq2 = AminoAcidSequence.([struc1[\"A\"], struc2[\"A\"]], standardselector, gaps=false)\n2-element Array{BioSequences.BioSequence{BioSequences.AminoAcidAlphabet},1}:\n MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLES…SHVKNLDENGLDLLSKMLIYDPAKRISGKMALNHPYFND\n MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRTEG…RSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL\n\njulia> using BioAlignments\n\njulia> scoremodel = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1);\n\njulia> alres = pairalign(GlobalAlignment(), seq1, seq2, scoremodel)\nPairwiseAlignmentResult{Int64,BioSequences.BioSequence{BioSequences.AminoAcidAlphabet},BioSequences.BioSequence{BioSequences.AminoAcidAlphabet}}:\n  score: 945\n  seq:   1 MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRH  60\n           ||   | ||||||||||||| | | || ||| |||| |    |||||||||||||||| |\n  ref:   1 MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRTE----GVPSTAIREISLLKELNH  56\n\n  seq:  61 PNIVSLQDVLMQDSRLYLIFEFLSMDLKKYLD-SIPPGQYMDSSLVKSYLYQILQGIVFC 119\n           |||| | ||      ||| ||||  ||||  | |   |      | |||| | |||  ||\n  ref:  57 PNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTG--IPLPLIKSYLFQLLQGLAFC 114\n\n  seq: 120 HSRRVLHRDLKPQNLLIDDKGTIKLADFGLARAFGV----YTHEVVTLWYRSPEVLLGSA 175\n           || ||||||||||||||   | ||||||||||||||    ||||||||||| || |||\n  ref: 115 HSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCK 174\n\n  seq: 176 RYSTPVDIWSIGTIFAELATKKPLFHGDSEIDQLFRIFRALGTPNNEVWPEVESLQDYKN 235\n            ||| ||||| | ||||  |   || ||||||||||||| ||||   ||| | |  |||\n  ref: 175 YYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKP 234\n\n  seq: 236 TFPKWKPGSLASHVKNLDENGLDLLSKMLIYDPAKRISGKMALNHPYFND---------- 285\n            ||||        |  ||| |  ||| || ||| |||| | || || | |\n  ref: 235 SFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL 294\nIn fact, pairalign is extended to carry out the above steps and return the alignment by calling pairalign(struc1[\"A\"], struc2[\"A\"], standardselector) in this case. scoremodel and aligntype are keyword arguments with the defaults shown above."
},

{
    "location": "documentation.html#Spatial-calculations-1",
    "page": "Documentation",
    "title": "Spatial calculations",
    "category": "section",
    "text": "Various functions are provided to calculate spatial quantities for proteins:Command Returns\ndistance Minimum distance between two elements\nsqdistance Minimum square distance between two elements\ncoordarray Atomic coordinates in Å of an element as a 2D Array with each column corresponding to one atom\nbondangle Angle between three atoms\ndihedralangle Dihedral angle defined by four atoms\nomegaangle Omega dihedral angle between a residue and the previous residue\nphiangle Phi dihedral angle between a residue and the previous residue\npsiangle Psi dihedral angle between a residue and the next residue\nomegaangles Vector of omega dihedral angles of an element\nphiangles Vector of phi dihedral angles of an element\npsiangles Vector of psi dihedral angles of an element\nramachandranangles Vectors of phi and psi angles of an element\nContactMap ContactMap of two elements, or one element with itself\nDistanceMap DistanceMap of two elements, or one element with itself\nshowcontactmap Print a representation of a ContactMap to stdout or a specified IO instance\nTransformation The 3D transformation to map one set of coordinates onto another\napplytransform! Modify all coordinates in an element according to a transformation\napplytransform Modify coordinates according to a transformation\nsuperimpose! Superimpose one element onto another\nrmsd RMSD between two elements, with or without superimposition\ndisplacements Vector of displacements between two elements, with or without superimposition\nMetaGraph Construct a MetaGraph of contacting elementsThe omegaangle, phiangle and psiangle functions can take either a pair of residues or a chain and a position. The omegaangle and phiangle functions measure the angle between the residue at the given index and the one before. The psiangle function measures between the given index and the one after. For example:julia> distance(struc[\'A\'][10], struc[\'A\'][20])\n10.782158874733762\n\njulia> rad2deg(bondangle(struc[\'A\'][50][\"N\"], struc[\'A\'][50][\"CA\"], struc[\'A\'][50][\"C\"]))\n110.77765846083398\n\njulia> rad2deg(dihedralangle(struc[\'A\'][50][\"N\"], struc[\'A\'][50][\"CA\"], struc[\'A\'][50][\"C\"], struc[\'A\'][51][\"N\"]))\n-177.38288114072924\n\njulia> rad2deg(psiangle(struc[\'A\'][50], struc[\'A\'][51]))\n-177.38288114072924\n\njulia> rad2deg(psiangle(struc[\'A\'], 50))\n-177.38288114072924ContactMap takes in a structural element or a list, such as a Chain or Vector{Atom}, and returns a ContactMap object showing the contacts between the elements for a specified distance. ContactMap can also be given two structural elements as arguments, in which case a non-symmetrical 2D array is returned showing contacts between the elements. The underlying BitArray{2} for ContactMap contacts can be accessed with contacts.data if required.julia> contacts = ContactMap(collectatoms(struc[\'A\'], cbetaselector), 8.0)\nContact map of size (85, 85)A plot recipe is defined for this so it can shown with Plots.jl:using Plots\nplot(contacts)(Image: contactmap)For a quick, text-based representation of a ContactMap use showcontactmap.DistanceMap works in an analogous way to ContactMap and gives a map of the distances. It can also be plotted:dists = DistanceMap(collectatoms(struc[\'A\'], cbetaselector))\nusing Plots\nplot(dists)(Image: distancemap)Structural elements can be superimposed, and superposition-dependent properties such as the RMSD can be calculated. To carry out superimposition, BioStructures.jl carries out a sequence alignment and superimposes aligned residues using the Kabsch algorithm. For example:# Change the coordinates of element 1 to superimpose it onto element 2\n# Do sequence alignment with standard residues and calculate the transformation with Cα atoms (the default)\nsuperimpose!(el1, el2, standardselector)\n\n# The transformation object for the above superimposition\nTransformation(el1, el2, standardselector)\n\n# Calculate the transformation with backbone atoms\nsuperimpose!(el1, el2, standardselector, alignatoms=backboneselector)\n\n# Calculate RMSD on Cα atoms (the default) after superimposition\nrmsd(el1, el2, standardselector)\n\n# Superimpose based on backbone atoms and calculate RMSD based on Cβ atoms\nrmsd(el1, el2, standardselector, alignatoms=backboneselector, rmsdatoms=cbetaselector)\n\n# Do not do a superimposition - assumes the elements are already superimposed\nrmsd(el1, el2, standardselector, superimpose=false)displacements is used in a similar way to rmsd but returns the vector of distances for each superimposed atom. dispatoms selects the atoms to calculate the displacements on.These transformation functions may be useful beyond the context of protein structures. For example, Transformation(c1, c2) calculates the transformation to map one set of coordinates to another. The coordinate sets must be the same size and have the number of dimensions in the first axis and the number of points in the second axis.The contacting elements in a molecular structure form a graph, and this can be retrieved using MetaGraph. This extends MetaGraph from MetaGraphs.jl, allowing you to use all the graph analysis tools in LightGraphs.jl. For example:julia> mg = MetaGraph(collectatoms(struc[\"A\"], cbetaselector), 8.0)\n{85, 423} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)\n\njulia> using LightGraphs, MetaGraphs\n\njulia> nv(mg)\n85\n\njulia> ne(mg)\n423\n\njulia> get_prop(mg, :contactdist)\n8.0\n\njulia> mg[10, :element]\nAtom CB with serial 71, coordinates [-3.766, 4.031, 23.526]See the LightGraphs docs for details on how to calculate properties such as shortest paths, centrality measures, community detection and more. Similar to ContactMap, contacts are found between any element type passed in. So if you wanted the graph of chain contacts in a protein complex you could give a Model as the first argument."
},

{
    "location": "documentation.html#Downloading-PDB-files-1",
    "page": "Documentation",
    "title": "Downloading PDB files",
    "category": "section",
    "text": "To download a PDB file to a specified directory:downloadpdb(\"1EN2\", dir=\"path/to/pdb/directory\")To download multiple PDB files to a specified directory:downloadpdb([\"1EN2\", \"1ALW\", \"1AKE\"], dir=\"path/to/pdb/directory\")To download a PDB file in PDB, XML, mmCIF or MMTF format use the format argument:# To get mmCIF\ndownloadpdb(\"1ALW\", dir=\"path/to/pdb/directory\", format=MMCIF)\n\n# To get XML\ndownloadpdb(\"1ALW\", dir=\"path/to/pdb/directory\", format=PDBXML)To apply a function to a downloaded file and delete the file afterwards:downloadpdb(f, \"1ALW\")Or, using Julia\'s do syntax:downloadpdb(\"1ALW\") do fp\n    s = read(fp, PDB)\n    # Do something\nendNote that some PDB entries, e.g. large viral assemblies, are not available as PDB format files. In this case consider downloading the mmCIF file or MMTF file instead."
},

{
    "location": "documentation.html#Reading-PDB-files-1",
    "page": "Documentation",
    "title": "Reading PDB files",
    "category": "section",
    "text": "To parse an existing PDB file into a Structure-Model-Chain-Residue-Atom framework:julia> struc = read(\"/path/to/pdb/file.pdb\", PDB)\nProteinStructure 1EN2.pdb with 1 models, 1 chains (A), 85 residues, 754 atomsRead a mmCIF/MMTF file instead by replacing PDB with MMCIF/MMTF. Various options can be set through optional keyword arguments when parsing PDB/mmCIF/MMTF files:Keyword Argument Description\nstructure_name::AbstractString The name given to the returned ProteinStructure; defaults to the file name\nremove_disorder::Bool=false Whether to remove atoms with alt loc ID not \' \' or \'A\'\nread_std_atoms::Bool=true Whether to read standard ATOM records\nread_het_atoms::Bool=true Whether to read HETATOM records\ngzip::Bool=false Whether the file is gzipped (MMTF files only)Use retrievepdb to download and parse a PDB file into a Structure-Model-Chain-Residue-Atom framework in a single line:julia> struc = retrievepdb(\"1ALW\", dir=\"path/to/pdb/directory\")\nINFO: Downloading PDB: 1ALW\nProteinStructure 1ALW.pdb with 1 models, 2 chains (A,B), 346 residues, 2928 atomsIf you prefer to work with data frames rather than the data structures in BioStructures, the DataFrame constructor from DataFrames.jl has been extended to construct relevant data frames from lists of atoms or residues:julia> df = DataFrame(collectatoms(struc));\n\njulia> first(df, 3)\n3×17 DataFrame. Omitted printing of 5 columns\n│ Row │ ishetero │ serial │ atomname │ altlocid │ resname │ chainid │ resnumber │ inscode │ x       │ y       │ z       │ occupancy │\n│     │ Bool     │ Int64  │ String   │ Char     │ String  │ String  │ Int64     │ Char    │ Float64 │ Float64 │ Float64 │ Float64   │\n├─────┼──────────┼────────┼──────────┼──────────┼─────────┼─────────┼───────────┼─────────┼─────────┼─────────┼─────────┼───────────┤\n│ 1   │ false    │ 1      │ N        │ \' \'      │ GLU     │ A       │ 94        │ \' \'     │ 15.637  │ -47.066 │ 18.179  │ 1.0       │\n│ 2   │ false    │ 2      │ CA       │ \' \'      │ GLU     │ A       │ 94        │ \' \'     │ 14.439  │ -47.978 │ 18.304  │ 1.0       │\n│ 3   │ false    │ 3      │ C        │ \' \'      │ GLU     │ A       │ 94        │ \' \'     │ 14.141  │ -48.183 │ 19.736  │ 1.0       │\n\njulia> df = DataFrame(collectresidues(struc));\n\njulia> first(df, 3)\n3×8 DataFrame\n│ Row │ ishetero │ resname │ chainid │ resnumber │ inscode │ countatoms │ modelnumber │ isdisorderedres │\n│     │ Bool     │ String  │ String  │ Int64     │ Char    │ Int64      │ Int64       │ Bool            │\n├─────┼──────────┼─────────┼─────────┼───────────┼─────────┼────────────┼─────────────┼─────────────────┤\n│ 1   │ false    │ GLU     │ A       │ 94        │ \' \'     │ 9          │ 1           │ false           │\n│ 2   │ false    │ GLU     │ A       │ 95        │ \' \'     │ 9          │ 1           │ false           │\n│ 3   │ false    │ VAL     │ A       │ 96        │ \' \'     │ 7          │ 1           │ false           │As with file writing disordered entities are expanded by default but this can be changed by setting expand_disordered to false."
},

{
    "location": "documentation.html#Writing-PDB-files-1",
    "page": "Documentation",
    "title": "Writing PDB files",
    "category": "section",
    "text": "PDB format files can be written:writepdb(\"1EN2_out.pdb\", struc)Any element type can be given as input to writepdb. The first argument can also be a stream. Atom selectors can also be given as additional arguments:# Only backbone atoms are written out\nwritepdb(\"1EN2_out.pdb\", struc, backboneselector)To write mmCIF format files, use the writemmcif function with similar arguments. A MMCIFDict can also be written using writemmcif:writemmcif(\"1EN2_out.dic\", mmcif_dict)To write out a MMTF file, use the writemmtf function with any element type or a MMTFDict as an argument. The gzip keyword argument, default false, determines whether to gzip the written file.Unlike for the collect functions, expand_disordered is set to true when writing files as it is usually desirable to retain all entities. Set expand_disordered to false to not write out more than one atom or residue at each location.Multi-character chain IDs can be written to mmCIF and MMTF files but will throw an error when written to a PDB file as the PDB file format only has one character allocated to the chain ID.If you want the PDB record line for an Atom, use pdbline. For example:julia> pdbline(at)\n\"HETATM  101  C  A  X B  20      10.500  20.123  -5.123  0.50 50.13           C1+\"If you want to generate a PDB record line from values directly, do so using an AtomRecord:julia> pdbline(AtomRecord(false, 669, \"CA\", \' \', \"ILE\", \"A\", 90, \' \', [31.743, 33.11, 31.221], 1.00, 25.76, \"C\", \"\"))\n\"ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  \"This can be useful when writing PDB files from your own data structures."
},

{
    "location": "documentation.html#RCSB-PDB-utility-functions-1",
    "page": "Documentation",
    "title": "RCSB PDB utility functions",
    "category": "section",
    "text": "To get the list of all PDB entries:l = pdbentrylist()To download the entire RCSB PDB database in your preferred file format:downloadentirepdb(dir=\"path/to/pdb/directory\", format=MMTF)This operation takes a lot of disk space and time to complete (depending on internet connection).To update your local PDB directory based on the weekly status list of new, modified and obsolete PDB files from the RCSB server:updatelocalpdb(dir=\"path/to/pdb/directory\", format=MMTF)Obsolete PDB files are stored in the auto-generated obsolete directory inside the specified local PDB directory.To maintain a local copy of the entire RCSB PDB database, run the downloadentirepdb function once to download all PDB files and set up a CRON job or similar to run updatelocalpdb function once a week to keep the local PDB directory up to date with the RCSB server.There are a few more functions that may be useful:Function Returns Return type\npdbstatuslist List of PDB entries from a specified RCSB weekly status list URL Array{String,1}\npdbrecentchanges Added, modified and obsolete PDB lists from the recent RCSB weekly status files Tuple{Array{String,1},Array{String,1}, Array{String,1}}\npdbobsoletelist List of all obsolete PDB entries Array{String,1}\ndownloadallobsoletepdb Downloads all obsolete PDB files from the RCSB PDB server Array{String,1}"
},

{
    "location": "documentation.html#Visualising-structures-1",
    "page": "Documentation",
    "title": "Visualising structures",
    "category": "section",
    "text": "The Bio3DView.jl package can be used to visualise molecular structures. For example:using Bio3DView\nusing Blink\nviewpdb(\"1CRN\")(Image: viewpdb)struc = retrievepdb(\"1AKE\")\nviewstruc(struc[\'A\'], surface=Surface(Dict(\"colorscheme\"=> \"greenCarbon\")))(Image: viewstruc)Here they are shown as static images but they are interactive when using Bio3DView.jl. See the Bio3DView.jl tutorial for more information."
},

{
    "location": "documentation.html#Related-software-1",
    "page": "Documentation",
    "title": "Related software",
    "category": "section",
    "text": "Other packages in the Julia ecosystem that deal with structural bioinformatics or related fields include:MIToS.jl - protein sequence and structure analysis.\nBio3DView.jl - view molecular structures (see Visualising structures).\nMMTF.jl - read and write MMTF files. BioStructures.jl builds on top of MMTF.jl.\nProteinEnsembles.jl - modelling ensembles of protein structures.\nMolly.jl - proof of concept molecular dynamics."
},

{
    "location": "examples.html#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples.html#BioStructures-examples-1",
    "page": "Examples",
    "title": "BioStructures examples",
    "category": "section",
    "text": "Further examples of BioStructures usage are given below.A) Plot the temperature factors of a protein:using Plots\ncalphas = collectatoms(struc, calphaselector)\nplot(resnumber.(calphas),\n     tempfactor.(calphas),\n     xlabel=\"Residue number\",\n     ylabel=\"Temperature factor\",\n     label=\"\")B) Print the PDB records for all Cα atoms within 5 Å of residue 38:for at in calphas\n    if distance(struc[\'A\'][38], at) < 5.0 && resnumber(at) != 38\n        println(pdbline(at))\n    end\nendC) Find the residues at the interface of a protein-protein interaction:for res_a in collectresidues(struc[\"A\"], standardselector)\n    for res_b in collectresidues(struc[\"B\"], standardselector)\n        if distance(res_a, res_b) < 5.0\n            println(resnumber(res_a), \"A \", resnumber(res_b), \"B\")\n        end\n    end\nendD) Show the Ramachandran phi/psi angle plot of a structure:using Plots\nphi_angles, psi_angles = ramachandranangles(struc, standardselector)\nscatter(rad2deg.(phi_angles),\n     rad2deg.(psi_angles),\n     title=\"Ramachandran plot\",\n     xlabel=\"Phi / degrees\",\n     ylabel=\"Psi / degrees\",\n     label=\"\",\n     xticks=[-180, -90, 0, 90, 180],\n     yticks=[-180, -90, 0, 90, 180],\n     xlims=(-180, 180),\n     ylims=(-180, 180))E) Calculate the RMSD and displacements between the heavy (non-hydrogen) atoms of two models in an NMR structure:downloadpdb(\"1SSU\")\nstruc_nmr = read(\"1SSU.pdb\", PDB)\nrmsd(struc_nmr[5], struc_nmr[10], superimpose=false, rmsdatoms=heavyatomselector)\ndisplacements(struc_nmr[5], struc_nmr[10], superimpose=false, rmsdatoms=heavyatomselector)F) Calculate the cysteine fraction of every structure in the PDB:l = pdbentrylist()\nfor p in l\n    downloadpdb(p, format=MMCIF) do fp\n        s = read(fp, MMCIF)\n        nres = countresidues(s, standardselector)\n        if nres > 0\n            frac = countresidues(s, standardselector, x -> resname(x) == \"CYS\") / nres\n            println(p, \"  \", round(frac, digits=2))\n        end\n    end\nendG) Interoperability is possible with other packages in the Julia ecosystem. For example, use NearestNeighbors.jl to find the 10 nearest residues to each residue:using NearestNeighbors\nstruc = retrievepdb(\"1AKE\")\nca = coordarray(struc[\"A\"], cbetaselector)\nkdtree = KDTree(ca; leafsize=10)\nidxs, dists = knn(kdtree, ca, 10, true)H) Interoperability with DataFrames.jl gives access to filtering, sorting, summary statistics and other writing options:using DataFrames\nusing CSV\nusing Statistics\nstruc = retrievepdb(\"1ALW\")\ndf = DataFrame(collectatoms(struc))\ndescribe(df) # Show summary\nmean(df.tempfactor) # Column-wise operations\nsort(df, :x) # Sorting\nCSV.write(\"1ALW.csv\", df) # CSV file writing"
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#BioStructures.StructuralElementOrList",
    "page": "API",
    "title": "BioStructures.StructuralElementOrList",
    "category": "constant",
    "text": "A StructuralElement or Vector of StructuralElements up to a Vector{Model}.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.backboneatomnames",
    "page": "API",
    "title": "BioStructures.backboneatomnames",
    "category": "constant",
    "text": "Set of protein backbone atom names.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.calphaatomnames",
    "page": "API",
    "title": "BioStructures.calphaatomnames",
    "category": "constant",
    "text": "Set of Cα atom names.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.cbetaatomnames",
    "page": "API",
    "title": "BioStructures.cbetaatomnames",
    "category": "constant",
    "text": "Set of Cβ atom names.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.pdbextension",
    "page": "API",
    "title": "BioStructures.pdbextension",
    "category": "constant",
    "text": "Mapping of Protein Data Bank (PDB) formats to their file extensions.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.waterresnames",
    "page": "API",
    "title": "BioStructures.waterresnames",
    "category": "constant",
    "text": "Set of residue names corresponding to water.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.AbstractAtom",
    "page": "API",
    "title": "BioStructures.AbstractAtom",
    "category": "type",
    "text": "An atom that is part of a macromolecule - either an Atom or a DisorderedAtom.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.AbstractResidue",
    "page": "API",
    "title": "BioStructures.AbstractResidue",
    "category": "type",
    "text": "A residue (amino acid) or other molecule - either a Residue or a DisorderedResidue.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.Atom",
    "page": "API",
    "title": "BioStructures.Atom",
    "category": "type",
    "text": "An atom that is part of a macromolecule.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.AtomRecord",
    "page": "API",
    "title": "BioStructures.AtomRecord",
    "category": "type",
    "text": "A record for a single atom, e.g. as represented in a Protein Data Bank (PDB) file.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.Chain",
    "page": "API",
    "title": "BioStructures.Chain",
    "category": "type",
    "text": "A chain (molecule) from a macromolecular structure.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.ContactMap",
    "page": "API",
    "title": "BioStructures.ContactMap",
    "category": "type",
    "text": "ContactMap(element, contact_distance)\nContactMap(element_one, element_two, contact_distance)\nContactMap(bit_array_2D)\n\nCalculate the contact map for a StructuralElementOrList, or between two StructuralElementOrLists.\n\nThis returns a ContactMap type containing a BitArray{2} with true where the sub-elements are no further than the contact distance and false otherwise. When one element is given as input this returns a symmetric square matrix.\n\nExamples\n\ncbetas_A = collectatoms(struc[\"A\"], cbetaselector)\ncbetas_B = collectatoms(struc[\"B\"], cbetaselector)\n\n# Contact map of chain A using standard C-beta and 8.0 Å definitions\nContactMap(cbetas_A, 8.0)\n\n# Rectangular contact map of chains A and B\nContactMap(cbetas_A, cbetas_B, 8.0)\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.DisorderedAtom",
    "page": "API",
    "title": "BioStructures.DisorderedAtom",
    "category": "type",
    "text": "A container to hold different locations of the same atom.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.DisorderedResidue",
    "page": "API",
    "title": "BioStructures.DisorderedResidue",
    "category": "type",
    "text": "A container to hold different versions of the same residue (point mutations).\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.DistanceMap",
    "page": "API",
    "title": "BioStructures.DistanceMap",
    "category": "type",
    "text": "DistanceMap(element)\nDistanceMap(element_one, element_two)\nDistanceMap(float_array_2D)\n\nCalculate the distance map for a StructuralElementOrList, or between two StructuralElementOrLists.\n\nThis returns a DistanceMap type containing a Array{Float64, 2} with minimum distances between the sub-elements. When one element is given as input this returns a symmetric square matrix.\n\nExamples\n\ncbetas_A = collectatoms(struc[\"A\"], cbetaselector)\ncbetas_B = collectatoms(struc[\"B\"], cbetaselector)\n\n# Distance map of chain A showing how far each residue is from the others\nDistanceMap(cbetas_A)\n\n# Rectangular distance map of chains A and B\nDistanceMap(cbetas_A, cbetas_B)\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.MMCIF",
    "page": "API",
    "title": "BioStructures.MMCIF",
    "category": "type",
    "text": "Protein Data Bank (PDB) mmCIF file format.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.MMCIFDict",
    "page": "API",
    "title": "BioStructures.MMCIFDict",
    "category": "type",
    "text": "MMCIFDict(filepath)\nMMCIFDict(io)\nMMCIFDict()\n\nA macromolecular Crystallographic Information File (mmCIF) dictionary.\n\nCan be accessed using similar functions to a standard Dict. Keys are field names as a String and values are always Vector{String}, even for multiple components or numerical data. To directly access the underlying dictionary of MMCIFDict d, use d.dict. Call MMCIFDict with a filepath or stream to read the dictionary from that source.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.MMTF",
    "page": "API",
    "title": "BioStructures.MMTF",
    "category": "type",
    "text": "Protein Data Bank (PDB) MMTF file format.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.MMTFDict",
    "page": "API",
    "title": "BioStructures.MMTFDict",
    "category": "type",
    "text": "MMTFDict(filepath)\nMMTFDict(io)\nMMTFDict()\n\nA Macromolecular Transmission Format (MMTF) dictionary.\n\nCan be accessed using similar functions to a standard Dict. Keys are field names as a String and values are various types. To directly access the underlying dictionary of MMTFDict d, use d.dict. Call MMTFDict with a filepath or stream to read the dictionary from that source. The keyword argument gzip (default false) determines if the file is gzipped.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.Model",
    "page": "API",
    "title": "BioStructures.Model",
    "category": "type",
    "text": "A conformation of a macromolecular structure.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.PDB",
    "page": "API",
    "title": "BioStructures.PDB",
    "category": "type",
    "text": "Protein Data Bank (PDB) file format.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.PDBParseError",
    "page": "API",
    "title": "BioStructures.PDBParseError",
    "category": "type",
    "text": "Error arising from parsing a Protein Data Bank (PDB) file.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.PDBXML",
    "page": "API",
    "title": "BioStructures.PDBXML",
    "category": "type",
    "text": "Protein Data Bank (PDB) XML file format.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.ProteinStructure",
    "page": "API",
    "title": "BioStructures.ProteinStructure",
    "category": "type",
    "text": "A container for multiple Models that represents a Protein Data Bank (PDB) entry.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.Residue",
    "page": "API",
    "title": "BioStructures.Residue",
    "category": "type",
    "text": "A residue (amino acid) or other molecule.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.SpatialMap",
    "page": "API",
    "title": "BioStructures.SpatialMap",
    "category": "type",
    "text": "A map of a structural property, e.g. a ContactMap or a DistanceMap.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.StructuralElement",
    "page": "API",
    "title": "BioStructures.StructuralElement",
    "category": "type",
    "text": "A macromolecular structural element.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.Transformation",
    "page": "API",
    "title": "BioStructures.Transformation",
    "category": "type",
    "text": "Transformation(el1, el2, residue_selectors...)\nTransformation(coords1, coords2)\nTransformation(trans1, trans2, rot)\n\nA 3D transformation to map one set of coordinates onto another. Found using the Kabsch algorithm. When called with structural elements, carries out a pairwise alignment and superimposes on atoms from aligned residues. In this case, keyword arguments for pairwise alignment can be given, see pairalign. The residue selectors determine which residues to do the pairwise alignment on. The keyword argument alignatoms is an atom selector that selects the atoms to calculate the superimposition on (default calphaselector). Can also be called with two sets of coordinates of the same size, with the number of dimensions in the first axis and the number of points in the second axis.\n\nThe returned Transformation object consists of the mean coordinates of the first set, the mean coordinates of the second set, the rotation to map the first centred set onto the second centred set, and the indices of the aligned residues in the first and second elements if relevant.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioCore.distance-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioCore.distance",
    "category": "method",
    "text": "distance(element_one, element_two, atom_selectors...)\n\nGet the minimum distance in Å between two StructuralElementOrLists.\n\nAdditional arguments are atom selector functions - only atoms that return true from the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.allselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.allselector",
    "category": "method",
    "text": "allselector(at)\nallselector(res)\n\nTrivial selector that returns true for any AbstractAtom or AbstractResidue. Use it to select all atoms or residues.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.altlocid-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.altlocid",
    "category": "method",
    "text": "altlocid(at)\n\nGet the alternative location ID of an AbstractAtom as a Char.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.altlocids-Tuple{DisorderedAtom}",
    "page": "API",
    "title": "BioStructures.altlocids",
    "category": "method",
    "text": "altlocids(dis_at)\n\nGet the list of alternative location IDs in an AbstractAtom as a Vector{Char}, sorted by atom serial.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.applyselectors!-Tuple{Array{#s21,1} where #s21<:StructuralElement,Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.applyselectors!",
    "category": "method",
    "text": "applyselectors!(els, selectors...)\n\nRemoves from a Vector of StructuralElements all elements that do not return true from all the selector functions.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.applyselectors-Tuple{Array{#s22,1} where #s22<:StructuralElement,Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.applyselectors",
    "category": "method",
    "text": "applyselectors(els, selectors...)\n\nReturns a copy of a Vector of StructuralElements with all elements that do not return true from all the selector functions removed.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.applytransform!-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Transformation}",
    "page": "API",
    "title": "BioStructures.applytransform!",
    "category": "method",
    "text": "applytransform!(el, transformation)\n\nModify all coordinates in an element according to a transformation.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.applytransform-Tuple{Array{#s182,2} where #s182<:Real,Transformation}",
    "page": "API",
    "title": "BioStructures.applytransform",
    "category": "method",
    "text": "applytransform(coords, transformation)\n\nModify coordinates according to a transformation.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.atomid-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.atomid",
    "category": "method",
    "text": "atomid(at)\n\nGet a descriptive atom ID for an AbstractAtom as a Tuple of the form (full residue ID, residue name, atom name).\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.atomname-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.atomname",
    "category": "method",
    "text": "atomname(at; strip=true)\n\nGet the atom name of an AbstractAtom as a String. strip determines whether surrounding whitespace is stripped.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.atomnames-Tuple{Residue}",
    "page": "API",
    "title": "BioStructures.atomnames",
    "category": "method",
    "text": "atomnames(res; strip=true)\n\nGet the sorted list of AbstractAtoms in an AbstractResidue. strip determines whether surrounding whitespace is stripped.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.atomnameselector-Tuple{AbstractAtom,Set{String}}",
    "page": "API",
    "title": "BioStructures.atomnameselector",
    "category": "method",
    "text": "atomnameselector(at, atom_names; strip=true)\n\nDetermines if an AbstractAtom has its atom name in the given Set or Vector. strip determines whether surrounding whitespace is stripped from the atom name before it is checked in the list.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.atoms-Tuple{Residue}",
    "page": "API",
    "title": "BioStructures.atoms",
    "category": "method",
    "text": "atoms(res)\n\nReturn the dictionary of AbstractAtoms in an AbstractResidue.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.backboneselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.backboneselector",
    "category": "method",
    "text": "backboneselector(at)\n\nDetermines if an AbstractAtom is not a hetero-atom and corresponds to a protein backbone atom.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.bondangle-Tuple{AbstractAtom,AbstractAtom,AbstractAtom}",
    "page": "API",
    "title": "BioStructures.bondangle",
    "category": "method",
    "text": "bondangle(atom_a, atom_b, atom_c)\nbondangle(vec_ba, vec_bc)\n\nCalculate the bond or pseudo-bond angle in radians between three AbstractAtoms or two vectors.\n\nThe angle between B→A and B→C is returned in the range 0 to π.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.calphaselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.calphaselector",
    "category": "method",
    "text": "calphaselector(at)\n\nDetermines if an AbstractAtom is not a hetero-atom and corresponds to a Cα atom.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.cbetaselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.cbetaselector",
    "category": "method",
    "text": "cbetaselector(at)\n\nDetermines if an AbstractAtom is not a hetero-atom and corresponds to a Cβ atom, or a Cα atom in glycine.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.chain-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.chain",
    "category": "method",
    "text": "chain(at)\nchain(res)\n\nReturn the Chain that an AbstractAtom or AbstractResidue belongs to.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.chainid-Tuple{Union{AbstractAtom, AbstractResidue}}",
    "page": "API",
    "title": "BioStructures.chainid",
    "category": "method",
    "text": "chainid(el)\n\nGet the chain ID of an AbstractAtom, AbstractResidue or Chain as a String.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.chainids-Tuple{Model}",
    "page": "API",
    "title": "BioStructures.chainids",
    "category": "method",
    "text": "chainids(mod)\nchainids(struc)\n\nGet the sorted chain IDs of the chains in a Model, or the default Model of a ProteinStructure, as a Vector{String}.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.chains-Tuple{Model}",
    "page": "API",
    "title": "BioStructures.chains",
    "category": "method",
    "text": "chains(mod)\nchains(struc)\n\nReturn the dictionary of Chains in a Model, or the default Model of a ProteinStructure.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.charge-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.charge",
    "category": "method",
    "text": "charge(at; strip=true)\n\nGet the charge on an AbstractAtom as a String. The charge is set to \"  \" if not specified during atom creation. strip determines whether surrounding whitespace is stripped.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.choosedefaultaltlocid-Tuple{Atom,Atom}",
    "page": "API",
    "title": "BioStructures.choosedefaultaltlocid",
    "category": "method",
    "text": "choosedefaultaltlocid(at_one, at_two)\n\nDetermine which of two Atoms representing a disorered atom better qualifies as the default location. The Atom with the highest occupancy is chosen; in the case of ties the Atom with the lowest alternative location ID in alphabetical order is chosen.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.collectatoms-Tuple{ProteinStructure}",
    "page": "API",
    "title": "BioStructures.collectatoms",
    "category": "method",
    "text": "collectatoms(el)\n\nReturns a Vector of the atoms in a StructuralElementOrList. Additional arguments are atom selector functions - only atoms that return true from all the functions are retained. The keyword argument expand_disordered (default false) determines whether to return all copies of disordered atoms separately.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.collectchains-Tuple{ProteinStructure}",
    "page": "API",
    "title": "BioStructures.collectchains",
    "category": "method",
    "text": "collectchains(el)\n\nReturns a Vector of the chains in a StructuralElementOrList. Additional arguments are chain selector functions - only chains that return true from all the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.collectmodels-Tuple{ProteinStructure}",
    "page": "API",
    "title": "BioStructures.collectmodels",
    "category": "method",
    "text": "collectmodels(el)\n\nReturns a Vector of the models in a StructuralElementOrList. Additional arguments are model selector functions - only models that return true from all the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.collectresidues-Tuple{ProteinStructure}",
    "page": "API",
    "title": "BioStructures.collectresidues",
    "category": "method",
    "text": "collectresidues(el)\n\nReturns a Vector of the residues in a StructuralElementOrList. Additional arguments are residue selector functions - only residues that return true from all the functions are retained. The keyword argument expand_disordered (default false) determines whether to return all copies of disordered residues separately.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.coordarray-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.coordarray",
    "category": "method",
    "text": "coordarray(element, atom_selectors...)\n\nGet the atomic coordinates in Å of a StructuralElementOrList as a 2D Array with each column corresponding to one atom.\n\nAdditional arguments are atom selector functions - only atoms that return true from all the functions are retained. The keyword argument expand_disordered (default false) determines whether to return coordinates for all copies of disordered atoms separately.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.coords!-Tuple{Atom,Array{#s22,1} where #s22<:Real}",
    "page": "API",
    "title": "BioStructures.coords!",
    "category": "method",
    "text": "coords!(at, new_coords)\n\nSet the coordinates of an AbstractAtom to a Vector of 3 numbers. For DisorderedAtoms only the default atom is updated.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.coords-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.coords",
    "category": "method",
    "text": "coords(at)\n\nGet the atomic coordinates of an AbstractAtom as a Vector{Float64}.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.countatoms-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.countatoms",
    "category": "method",
    "text": "countatoms(el)\n\nReturn the number of atoms in a StructuralElementOrList as an Int. Additional arguments are atom selector functions - only atoms that return true from all the functions are counted. The keyword argument expand_disordered (default false) determines whether to return all copies of disordered atoms separately.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.countchains-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.countchains",
    "category": "method",
    "text": "countchains(el)\n\nReturn the number of Chains in a StructuralElementOrList as an Int. Additional arguments are chain selector functions - only chains that return true from all the functions are counted.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.countmodels-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.countmodels",
    "category": "method",
    "text": "countmodels(el)\n\nReturn the number of Models in a StructuralElementOrList as an Int. Additional arguments are model selector functions - only models that return true from all the functions are counted.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.countresidues-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.countresidues",
    "category": "method",
    "text": "countresidues(el)\n\nReturn the number of residues in a StructuralElementOrList as an Int. Additional arguments are residue selector functions - only residues that return true from all the functions are counted. The keyword argument expand_disordered (default false) determines whether to return all copies of disordered residues separately.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.defaultaltlocid-Tuple{DisorderedAtom}",
    "page": "API",
    "title": "BioStructures.defaultaltlocid",
    "category": "method",
    "text": "defaultaltlocid(dis_at)\n\nGet the alternative location ID of the default Atom in a DisorderedAtom as a Char. The default is the highest occupancy, or lowest character alternative location ID for ties (i.e. \'A\' beats \'B\').\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.defaultatom-Tuple{DisorderedAtom}",
    "page": "API",
    "title": "BioStructures.defaultatom",
    "category": "method",
    "text": "defaultatom(dis_at)\n\nReturn the default Atom in a DisorderedAtom. The default is the highest occupancy, or lowest character alternative location ID for ties (i.e. \'A\' beats \'B\').\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.defaultmodel-Tuple{ProteinStructure}",
    "page": "API",
    "title": "BioStructures.defaultmodel",
    "category": "method",
    "text": "defaultmodel(struc)\n\nGet the default Model in a ProteinStructure. This is the Model with the lowest model number.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.defaultresidue-Tuple{DisorderedResidue}",
    "page": "API",
    "title": "BioStructures.defaultresidue",
    "category": "method",
    "text": "defaultresidue(dis_res)\n\nReturn the default Residue in a DisorderedResidue. The default is the first name read in.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.defaultresname-Tuple{DisorderedResidue}",
    "page": "API",
    "title": "BioStructures.defaultresname",
    "category": "method",
    "text": "defaultresname(dis_res)\n\nGet the name of the default Residue in a DisorderedResidue as a String. The default is the first name read in.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.dihedralangle-NTuple{4,AbstractAtom}",
    "page": "API",
    "title": "BioStructures.dihedralangle",
    "category": "method",
    "text": "dihedralangle(atom_a, atom_b, atom_c, atom_d)\ndihedralangle(vec_ab, vec_bc, vec_cd)\n\nCalculate the dihedral angle in radians defined by four AbstractAtoms or three vectors.\n\nThe angle between the planes defined by atoms (A, B, C) and (B, C, D) is returned in the range -π to π.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.disorderedres-Tuple{DisorderedResidue,AbstractString}",
    "page": "API",
    "title": "BioStructures.disorderedres",
    "category": "method",
    "text": "disorderedres(dis_res, res_name)\n\nReturn the Residue in a DisorderedResidue with a given residue name.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.disorderselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.disorderselector",
    "category": "method",
    "text": "disorderselector(at)\ndisorderselector(res)\n\nDetermines whether an AbstractAtom or AbstractResidue is disordered, i.e. has multiple locations in the case of atoms or multiple residue names (point mutants) in the case of residues.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.displacements-Tuple{Array{#s181,N} where N where #s181<:Real,Array{#s180,N} where N where #s180<:Real}",
    "page": "API",
    "title": "BioStructures.displacements",
    "category": "method",
    "text": "displacements(element_one, element_two, residue_selectors...)\ndisplacements(element_one, element_two, superimpose=false)\ndisplacements(coords_one, coords_two)\n\nGet the displacements in Å between atomic coordinates from two StructuralElementOrLists or two coordinate Arrays.\n\nIf superimpose is true (the default), the elements are superimposed before calculation and the displacements are calculated on the superimposed residues. See Transformation for keyword arguments. If superimpose is false the elements are assumed to be superimposed and must be of the same length. The keyword argument dispatoms is an atom selector that selects the atoms to calculate displacements on (default calphaselector).\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.downloadallobsoletepdb-Tuple{}",
    "page": "API",
    "title": "BioStructures.downloadallobsoletepdb",
    "category": "method",
    "text": "downloadallobsoletepdb(; kwargs...)\n\nDownload all obsolete Protein Data Bank (PDB) files from the RCSB server.\n\nReturns the list of PDB IDs downloaded. Requires an internet connection.\n\nArguments\n\nobsolete_dir::AbstractString=pwd(): the directory where the PDB files are   downloaded; defaults to the current working directory.\nformat::Type=PDB: the format of the PDB file; options are PDB, PDBXML,   MMCIF and MMTF.\noverwrite::Bool=false: if set true, overwrites the PDB file if it exists   in dir; by default skips downloading the PDB file if it exists.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.downloadentirepdb-Tuple{}",
    "page": "API",
    "title": "BioStructures.downloadentirepdb",
    "category": "method",
    "text": "downloadentirepdb(; kwargs...)\n\nDownload the entire Protein Data Bank (PDB) from the RCSB server.\n\nReturns the list of PDB IDs downloaded. Ensure you have enough disk space and time before running. The function can be stopped any time and called again to resume downloading. Requires an internet connection.\n\nArguments\n\ndir::AbstractString=pwd(): the directory to which the PDB files are   downloaded; defaults to the current working directory.\nformat::Type=PDB: the format of the PDB file; options are PDB, PDBXML,   MMCIF and MMTF.\noverwrite::Bool=false: if set true, overwrites the PDB file if it exists   in dir; by default skips downloading the PDB file if it exists.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.downloadpdb-Tuple{AbstractString}",
    "page": "API",
    "title": "BioStructures.downloadpdb",
    "category": "method",
    "text": "downloadpdb(pdbid::AbstractString; kwargs...)\ndownloadpdb(pdbid::Array{<:AbstractString, 1}; kwargs...)\ndownloadpdb(f::Function, args...)\n\nDownload files from the Protein Data Bank (PDB) via RCSB.\n\nWhen given an AbstractString, e.g. \"1AKE\", downloads the PDB file and returns the path to the file. When given an Array{<:AbstractString, 1}, downloads the PDB files in the array and returns an array of the paths to the files. When given a function as the first argument, runs the function with the downloaded filepath(s) as an argument then removes the file(s). Requires an internet connection.\n\nArguments\n\ndir::AbstractString=pwd(): the directory to which the PDB file is   downloaded; defaults to the current working directory.\nformat::Type=PDB: the format of the PDB file; options are PDB, PDBXML,   MMCIF and MMTF.\nobsolete::Bool=false: if set true, the PDB file is downloaded in the   auto-generated \"obsolete\" directory inside the specified dir.\noverwrite::Bool=false: if set true, overwrites the PDB file if it exists   in dir; by default skips downloading the PDB file if it exists.\nba_number::Integer=0: if set > 0 downloads the respective biological   assembly; by default downloads the PDB file.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.element-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.element",
    "category": "method",
    "text": "element(at; strip=true)\n\nGet the element of an AbstractAtom as a String. The element is set to \"  \" if not specified during atom creation. strip determines whether surrounding whitespace is stripped.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.generatechainid-Tuple{Integer}",
    "page": "API",
    "title": "BioStructures.generatechainid",
    "category": "method",
    "text": "generatechainid(entity_id)\n\nConvert a positive Integer into a chain ID. Goes A to Z, then AA to ZA, AB to ZB etc. This is in line with PDB conventions.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.heavyatomselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.heavyatomselector",
    "category": "method",
    "text": "heavyatomselector(at)\n\nDetermines if an AbstractAtom corresponds to a heavy (non-hydrogen) atom and is not a hetero-atom.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.heteroselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.heteroselector",
    "category": "method",
    "text": "heteroselector(at)\nheteroselector(res)\n\nDetermines if an AbstractAtom represents a hetero atom, e.g. came from a HETATM record in a Protein Data Bank (PDB) file, or if an AbstractResidue represents a hetero molecule, e.g. consists of HETATM records from a PDB file.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.hydrogenselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.hydrogenselector",
    "category": "method",
    "text": "hydrogenselector(at)\n\nDetermines if an AbstractAtom represents hydrogen. Uses the element field where possible, otherwise uses the atom name.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.inscode-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.inscode",
    "category": "method",
    "text": "inscode(at)\ninscode(res)\n\nGet the insertion code of an AbstractAtom or AbstractResidue as a Char.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.isdisorderedatom-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.isdisorderedatom",
    "category": "method",
    "text": "isdisorderedatom(at)\n\nDetermines if an AbstractAtom is a DisorderedAtom, i.e. if there are multiple locations present for an atom.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.isdisorderedres-Tuple{Residue}",
    "page": "API",
    "title": "BioStructures.isdisorderedres",
    "category": "method",
    "text": "isdisorderedres(res)\n\nDetermine if an AbstractResidue is a DisorderedResidue, i.e. there are multiple residue names with the same residue ID.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.ishetero-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.ishetero",
    "category": "method",
    "text": "ishetero(at)\nishetero(res)\n\nDetermines if an AbstractAtom represents a hetero atom, e.g. came from a HETATM record in a Protein Data Bank (PDB) file, or if an AbstractResidue represents a hetero molecule, e.g. consists of HETATM records from a PDB file.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.model-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.model",
    "category": "method",
    "text": "model(el)\n\nReturn the Model that an AbstractAtom, AbstractResidue or Chain belongs to.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.modelnumber-Tuple{Model}",
    "page": "API",
    "title": "BioStructures.modelnumber",
    "category": "method",
    "text": "modelnumber(el)\n\nGet the model number of a Model, Chain, AbstractResidue or AbstractAtom as an Int.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.modelnumbers-Tuple{ProteinStructure}",
    "page": "API",
    "title": "BioStructures.modelnumbers",
    "category": "method",
    "text": "modelnumbers(struc)\n\nGet the sorted model numbers from a ProteinStructure as a Vector{Int}.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.models-Tuple{ProteinStructure}",
    "page": "API",
    "title": "BioStructures.models",
    "category": "method",
    "text": "models(struc)\n\nReturn the dictionary of Models in a ProteinStructure.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.notwaterselector-Tuple{Union{AbstractAtom, AbstractResidue}}",
    "page": "API",
    "title": "BioStructures.notwaterselector",
    "category": "method",
    "text": "notwaterselector(res)\nnotwaterselector(at)\n\nDetermines if an AbstractResidue or AbstractAtom does not represent a water molecule.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.occupancy-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.occupancy",
    "category": "method",
    "text": "occupancy(at)\n\nGet the occupancy of an AbstractAtom as a Float64. The occupancy is set to 1.0 if not specified during atom creation.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.omegaangle-Tuple{AbstractResidue,AbstractResidue}",
    "page": "API",
    "title": "BioStructures.omegaangle",
    "category": "method",
    "text": "omegaangle(res, res_previous)\nomegaangle(chain, res_id)\n\nCalculate the omega angle in radians for an AbstractResidue.\n\nArguments can either be a residue and the previous residue or a chain and the position as a residue ID. The first residue (or one at the given index) requires the atoms \"N\" and \"CA\" and the previous residue requires the atoms \"CA\" and \"C\". The angle is in the range -π to π.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.omegaangles-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.omegaangles",
    "category": "method",
    "text": "omegaangles(element, residue_selectors...)\n\nCalculate the Vector of omega angles of a StructuralElementOrList.\n\nThe vectors have NaN for residues where an angle cannot be calculated, e.g. due to missing atoms or lack of an adjacent residue. The angle is in the range -π to π. Additional arguments are residue selector functions - only residues that return true from the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.pdbentrylist-Tuple{}",
    "page": "API",
    "title": "BioStructures.pdbentrylist",
    "category": "method",
    "text": "pdbentrylist()\n\nObtain the list of all Protein Data Bank (PDB) entries from the RCSB server.\n\nRequires an internet connection.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.pdbline-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.pdbline",
    "category": "method",
    "text": "pdbline(at::Atom)\npdbline(at::DisorderedAtom)\npdbline(at::AtomRecord)\n\nForm a Protein Data Bank (PDB) format ATOM or HETATM record as a String from an Atom, DisorderedAtom or AtomRecord.\n\nThis will throw an ArgumentError if a value cannot fit into the allocated space, e.g. the chain ID is longer than one character or the atom serial is greater than 99999. In this case consider using writemmcif or writemmtf to write a mmCIF file or a MMTF file.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.pdbobsoletelist-Tuple{}",
    "page": "API",
    "title": "BioStructures.pdbobsoletelist",
    "category": "method",
    "text": "pdbobsoletelist()\n\nObtain the list of all obsolete Protein Data Bank (PDB) entries from the RCSB server.\n\nRequires an internet connection.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.pdbrecentchanges-Tuple{}",
    "page": "API",
    "title": "BioStructures.pdbrecentchanges",
    "category": "method",
    "text": "pdbrecentchanges()\n\nObtain three lists giving the added, modified and obsolete Protein Data Bank (PDB) entries from the recent RCSB weekly status files.\n\nRequires an internet connection.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.pdbstatuslist-Tuple{AbstractString}",
    "page": "API",
    "title": "BioStructures.pdbstatuslist",
    "category": "method",
    "text": "pdbstatuslist(url::AbstractString)\n\nObtain the list of Protein Data Bank (PDB) entries from a RCSB weekly status file by specifying its URL.\n\nAn example URL is ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/added.pdb. Requires an internet connection.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.phiangle-Tuple{AbstractResidue,AbstractResidue}",
    "page": "API",
    "title": "BioStructures.phiangle",
    "category": "method",
    "text": "phiangle(res, res_previous)\nphiangle(chain, res_id)\n\nCalculate the phi angle in radians for an AbstractResidue.\n\nArguments can either be a residue and the previous residue or a chain and the position as a residue ID. The first residue (or one at the given index) requires the atoms \"N\", \"CA\" and \"C\" and the previous residue requires the atom \"C\". The angle is in the range -π to π.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.phiangles-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.phiangles",
    "category": "method",
    "text": "phiangles(element, residue_selectors...)\n\nCalculate the Vector of phi angles of a StructuralElementOrList.\n\nThe vectors have NaN for residues where an angle cannot be calculated, e.g. due to missing atoms or lack of an adjacent residue. The angle is in the range -π to π. Additional arguments are residue selector functions - only residues that return true from the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.psiangle-Tuple{AbstractResidue,AbstractResidue}",
    "page": "API",
    "title": "BioStructures.psiangle",
    "category": "method",
    "text": "psiangle(res, res_next)\npsiangle(chain, res_id)\n\nCalculate the psi angle in radians for an AbstractResidue.\n\nArguments can either be a residue and the next residue or a chain and the position as a residue ID. The first residue (or one at the given index) requires the atoms \"N\", \"CA\" and \"C\" and the next residue requires the atom \"N\". The angle is in the range -π to π.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.psiangles-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.psiangles",
    "category": "method",
    "text": "psiangles(element, residue_selectors...)\n\nCalculate the Vector of psi angles of a StructuralElementOrList.\n\nThe vectors have NaN for residues where an angle cannot be calculated, e.g. due to missing atoms or lack of an adjacent residue. The angle is in the range -π to π. Additional arguments are residue selector functions - only residues that return true from the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.ramachandranangles-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.ramachandranangles",
    "category": "method",
    "text": "ramachandranangles(element, residue_selectors...)\n\nCalculate the Vectors of phi and psi angles of a StructuralElementOrList.\n\nThe vectors have NaN for residues where an angle cannot be calculated, e.g. due to missing atoms or lack of an adjacent residue. The angles are in the range -π to π. Additional arguments are residue selector functions - only residues that return true from the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.resid-Tuple{AbstractResidue}",
    "page": "API",
    "title": "BioStructures.resid",
    "category": "method",
    "text": "resid(res; full=true)\n\nGet a descriptive residue ID String for an AbstractAtom or AbstractResidue. Format is residue number then insertion code with \"H\" in front for hetero residues. If full equals true the chain ID is also added after a colon. Examples are \"50A\", \"H20\" and \"10:A\".\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.resids-Tuple{Chain}",
    "page": "API",
    "title": "BioStructures.resids",
    "category": "method",
    "text": "resids(ch)\n\nGet the sorted list of AbstractResidues in a Chain.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.residue-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.residue",
    "category": "method",
    "text": "residue(at)\n\nGet the Residue that an AbstractAtom belongs to.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.residues-Tuple{Chain}",
    "page": "API",
    "title": "BioStructures.residues",
    "category": "method",
    "text": "residues(ch)\n\nReturn the dictionary of AbstractResidues in a Chain.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.resname-Tuple{Residue}",
    "page": "API",
    "title": "BioStructures.resname",
    "category": "method",
    "text": "resname(at; strip=true)\nresname(res; strip=true)\n\nGet the residue name of an AbstractAtom or AbstractResidue as a String. strip determines whether surrounding whitespace is stripped.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.resnames-Tuple{DisorderedResidue}",
    "page": "API",
    "title": "BioStructures.resnames",
    "category": "method",
    "text": "resnames(dis_res)\n\nGet the residue names in an AbstractResidue as a Vector{String}. For a DisorderedResidue there will be multiple residue names - in this case the default residue name is placed first, then the others are ordered alphabetically.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.resnameselector-Tuple{Union{AbstractAtom, AbstractResidue},Set{String}}",
    "page": "API",
    "title": "BioStructures.resnameselector",
    "category": "method",
    "text": "resnameselector(res, res_names)\nresnameselector(at, res_names)\n\nDetermines if an AbstractResidue or AbstractAtom has its residue name in the given Set or Vector.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.resnumber-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.resnumber",
    "category": "method",
    "text": "resnumber(at)\nresnumber(res)\n\nGet the residue number of an AbstractAtom or AbstractResidue as an Int.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.retrievepdb-Tuple{AbstractString}",
    "page": "API",
    "title": "BioStructures.retrievepdb",
    "category": "method",
    "text": "retrievepdb(pdbid::AbstractString; kwargs...)\n\nDownload and read a Protein Data Bank (PDB) file or biological assembly from the RCSB server, returning a ProteinStructure.\n\nRequires an internet connection.\n\nArguments\n\npdbid::AbstractString: the PDB ID to be downloaded and read.\ndir::AbstractString=pwd(): the directory to which the PDB file is   downloaded; defaults to the current working directory.\nobsolete::Bool=false: if set true, the PDB file is downloaded in the   auto-generated \"obsolete\" directory inside the specified dir.\noverwrite::Bool=false: if set true, overwrites the PDB file if it exists   in dir; by default skips downloading the PDB file if it exists.\nba_number::Integer=0: if set > 0 downloads the respective biological   assembly; by default downloads the PDB file.\nstructure_name::AbstractString=\"$pdbid.pdb\": the name given to the returned   ProteinStructure; defaults to the PDB ID.\nremove_disorder::Bool=false: whether to remove atoms with alt loc ID not \' \'   or \'A\'.\nread_std_atoms::Bool=true: whether to read standard ATOM records.\nread_het_atoms::Bool=true: whether to read HETATOM records.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.rmsd-Tuple{Array{#s181,N} where N where #s181<:Real,Array{#s180,N} where N where #s180<:Real}",
    "page": "API",
    "title": "BioStructures.rmsd",
    "category": "method",
    "text": "rmsd(element_one, element_two, residue_selectors...)\nrmsd(element_one, element_two, superimpose=false)\nrmsd(coords_one, coords_two)\n\nGet the root-mean-square deviation (RMSD) in Å between two StructuralElementOrLists or two coordinate Arrays.\n\nIf superimpose is true (the default), the elements are superimposed before RMSD calculation and the RMSD is calculated on the superimposed residues. See Transformation for keyword arguments. If superimpose is false the elements are assumed to be superimposed and must be of the same length. The keyword argument rmsdatoms is an atom selector that selects the atoms to calculate RMSD on (default calphaselector).\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.sequentialresidues-Tuple{AbstractResidue,AbstractResidue}",
    "page": "API",
    "title": "BioStructures.sequentialresidues",
    "category": "method",
    "text": "sequentialresidues(res_first, res_second)\n\nDetermine if the second residue follows the first in sequence. For this to be true the residues need to have the same chain ID, both need to be standard/hetero residues and the residue number of the second needs to be one greater than that of the first (or the residue numbers the same and the insertion code of the second greater than the first).\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.serial-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.serial",
    "category": "method",
    "text": "serial(at)\n\nGet the serial number of an AbstractAtom as an Int.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.showcontactmap-Tuple{IO,ContactMap}",
    "page": "API",
    "title": "BioStructures.showcontactmap",
    "category": "method",
    "text": "showcontactmap(contact_map)\nshowcontactmap(io, contact_map)\n\nPrint a representation of a ContactMap to stdout, or a specified IO instance. A fully plotted version can be obtained with plot(contact_map) but that requires Plots.jl; showcontactmap works without that dependency.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.spaceatomname-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.spaceatomname",
    "category": "method",
    "text": "spaceatomname(at::Atom)\n\nSpace an Atom name such that the last element letter (generally) appears in the second column.\n\nIf the element property of the Atom is set it is used to get the element, otherwise the name starts from the second column where possible. This function is generally not required as spacing is recorded when atom names are read in from a Protein Data Bank (PDB) file. However this spacing can be important, for example distinguising between Cα and calcium atoms.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.sqdistance-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.sqdistance",
    "category": "method",
    "text": "sqdistance(element_one, element_two, atom_selectors...)\n\nGet the minimum square distance in Å between two StructuralElementOrLists.\n\nAdditional arguments are atom selector functions - only atoms that return true from the functions are retained.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.standardselector-Tuple{AbstractAtom}",
    "page": "API",
    "title": "BioStructures.standardselector",
    "category": "method",
    "text": "standardselector(at)\nstandardselector(res)\n\nDetermines if an AbstractAtom represents a standard atom, e.g. came from a ATOM record in a Protein Data Bank (PDB) file, or if an AbstractResidue represents a standard molecule, e.g. consists of ATOM records from a PDB file.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.structure-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.structure",
    "category": "method",
    "text": "structure(el)\n\nReturn the ProteinStructure that an AbstractAtom, AbstractResidue, Chain or Model belongs to.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.structurename-Tuple{Union{AbstractAtom, AbstractResidue, Chain, Model}}",
    "page": "API",
    "title": "BioStructures.structurename",
    "category": "method",
    "text": "structurename(el)\n\nGet the name of the ProteinStructure that a StructuralElement belongs to as a String.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.superimpose!-Tuple{Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.superimpose!",
    "category": "method",
    "text": "superimpose!(el1, el2, residue_selectors...)\n\nCalculate the Transformation that maps the first element onto the second, and modify all coordinates in the first element according to the transformation. See Transformation for keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.tempfactor-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.tempfactor",
    "category": "method",
    "text": "tempfactor(at)\n\nGet the temperature factor of an AbstractAtom as a Float64. The temperature factor is set to 0.0 if not specified during atom creation.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.updatelocalpdb-Tuple{}",
    "page": "API",
    "title": "BioStructures.updatelocalpdb",
    "category": "method",
    "text": "updatelocalpdb(; dir::AbstractString=pwd(), format::Type=PDB)\n\nUpdate a local copy of the Protein Data Bank (PDB).\n\nObtains the recent weekly lists of new, modified and obsolete PDB entries and automatically updates the PDB files of the given format inside the local dir directory. Requires an internet connection.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.waterselector-Tuple{Union{AbstractAtom, AbstractResidue}}",
    "page": "API",
    "title": "BioStructures.waterselector",
    "category": "method",
    "text": "waterselector(res)\nwaterselector(at)\n\nDetermines if an AbstractResidue or AbstractAtom represents a water molecule.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.writemmcif-Tuple{AbstractString,MMCIFDict}",
    "page": "API",
    "title": "BioStructures.writemmcif",
    "category": "method",
    "text": "writemmcif(output, element, atom_selectors...)\nwritemmcif(output, mmcif_dict)\n\nWrite a StructuralElementOrList or a MMCIFDict to a mmCIF format file or output stream.\n\nAtom selector functions can be given as additional arguments - only atoms that return true from all the functions are retained. The keyword argument expand_disordered (default true) determines whether to return all copies of disordered residues and atoms.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.writepdb-Tuple{AbstractString,Union{StructuralElement, Array{Chain,1}, Array{Model,1}, Array{#s23,1} where #s23<:AbstractResidue, Array{#s22,1} where #s22<:AbstractAtom},Vararg{Function,N} where N}",
    "page": "API",
    "title": "BioStructures.writepdb",
    "category": "method",
    "text": "writepdb(output, element, atom_selectors...)\n\nWrite a StructuralElementOrList to a Protein Data Bank (PDB) format file or output stream.\n\nOnly ATOM, HETATM, MODEL and ENDMDL records are written - there is no header and there are no TER records. Atom selector functions can be given as additional arguments - only atoms that return true from all the functions are retained. The keyword argument expand_disordered (default true) determines whether to return all copies of disordered residues and atoms.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.x!-Tuple{Atom,Real}",
    "page": "API",
    "title": "BioStructures.x!",
    "category": "method",
    "text": "x!(at, val)\n\nSet the x coordinate of an AbstractAtom to val. For DisorderedAtoms only the default atom is updated.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.x-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.x",
    "category": "method",
    "text": "x(at)\n\nGet the x coordinate of an AbstractAtom as a Float64.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.y!-Tuple{Atom,Real}",
    "page": "API",
    "title": "BioStructures.y!",
    "category": "method",
    "text": "y!(at, val)\n\nSet the y coordinate of an AbstractAtom to val. For DisorderedAtoms only the default atom is updated.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.y-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.y",
    "category": "method",
    "text": "y(at)\n\nGet the y coordinate of an AbstractAtom as a Float64.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.z!-Tuple{Atom,Real}",
    "page": "API",
    "title": "BioStructures.z!",
    "category": "method",
    "text": "z!(at, val)\n\nSet the z coordinate of an AbstractAtom to val. For DisorderedAtoms only the default atom is updated.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures.z-Tuple{Atom}",
    "page": "API",
    "title": "BioStructures.z",
    "category": "method",
    "text": "z(at)\n\nGet the z coordinate of an AbstractAtom as a Float64.\n\n\n\n\n\n"
},

{
    "location": "api.html#MMTF.writemmtf-Tuple{Union{AbstractString, IO},MMTFDict}",
    "page": "API",
    "title": "MMTF.writemmtf",
    "category": "method",
    "text": "writemmtf(output, element, atom_selectors...)\nwritemmtf(output, mmtf_dict)\n\nWrite a StructuralElementOrList or a MMTFDict to a MMTF file or output stream.\n\nAtom selector functions can be given as additional arguments - only atoms that return true from all the functions are retained. The keyword argument expand_disordered (default true) determines whether to return all copies of disordered residues and atoms. The keyword argument gzip (default false) determines if the file should be gzipped.\n\n\n\n\n\n"
},

{
    "location": "api.html#BioStructures-API-1",
    "page": "API",
    "title": "BioStructures API",
    "category": "section",
    "text": "Modules = [BioStructures]\nPrivate = false"
},

]}
