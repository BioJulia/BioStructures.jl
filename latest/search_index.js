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
    "text": "BioStructures provides functionality to read, write and manipulate macromolecular structures, in particular proteins. Protein Data Bank (PDB) and mmCIF format files can be read in to a hierarchical data structure. Spatial calculations and functions to access the PDB are also provided."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install BioStructures from the Julia package REPL, which can be accessed by pressing ]:add BioStructuresIf you are interested in the cutting edge of the development, please check out the master branch to try new features before release."
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
    "text": "The BioStructures.jl package provides functionality to manipulate macromolecular structures, and in particular to read and write Protein Data Bank (PDB) and mmCIF files. It is designed to be used for standard structural analysis tasks, as well as acting as a platform on which others can build to create more specific tools.It compares favourably in terms of performance to other PDB parsers - see some benchmarks. The PDB and mmCIF parsers currently read in the whole PDB without explicit errors (with one exception). Help can be found on individual functions using ?function_name."
},

{
    "location": "documentation.html#Basics-1",
    "page": "Documentation",
    "title": "Basics",
    "category": "section",
    "text": "To download a PDB file:using BioStructures\n\n# Stored in the current working directory by default\ndownloadpdb(\"1EN2\")To parse a PDB file into a Structure-Model-Chain-Residue-Atom framework:julia> struc = read(\"/path/to/pdb/file.pdb\", PDB)\nProteinStructure 1EN2.pdb with 1 models, 1 chains (A), 85 residues, 754 atomsmmCIF files can be read into the same data structure with read(\"/path/to/cif/file.cif\", MMCIF). If you want to read an mmCIF file into a dictionary to query yourself (e.g. to access metadata fields), use MMCIFDict:julia> mmcif_dict = MMCIFDict(\"/path/to/cif/file.cif\")\nmmCIF dictionary with 716 fields\n\njulia> mmcif_dict[\"_entity_src_nat.common_name\"]\n1-element Array{String,1}:\n \"great nettle\"A MMCIFDict can be accessed in similar ways to a standard dictionary, and if necessary the underlying dictionary of MMCIFDict d can be accessed with d.dict. Note that the elements of the dictionary are always an Array{String,1}, even if only one value was read in or the data is numerical.Refer to Downloading PDB files and Reading PDB files sections for more options.The elements of struc can be accessed as follows:Command Returns Return type\nstruc[1] Model 1 Model\nstruc[1][\"A\"] Model 1, chain A Chain\nstruc[1][\'A\'] Shortcut to above if the chain ID is a single character Chain\nstruc[\"A\"] The lowest model (model 1), chain A Chain\nstruc[\"A\"][\"50\"] Model 1, chain A, residue 50 AbstractResidue\nstruc[\"A\"][50] Shortcut to above if it is not a hetero residue and the insertion code is blank AbstractResidue\nstruc[\"A\"][\"H_90\"] Model 1, chain A, hetero residue 90 AbstractResidue\nstruc[\"A\"][50][\"CA\"] Model 1, chain A, residue 50, atom name CA AbstractAtom\nstruc[\"A\"][15][\"CG\"][\'A\'] For disordered atoms, access a specific location AtomDisordered atoms are stored in a DisorderedAtom container but calls fall back to the default atom, so disorder can be ignored if you are not interested in it. Disordered residues (i.e. point mutations with different residue names) are stored in a DisorderedResidue container.The idea is that disorder will only bother you if you want it to. See the Biopython discussion for more.Properties can be retrieved as follows:Function Returns Return type\nserial Serial number of an atom Int\natomname Name of an atom String\naltlocid Alternative location ID of an atom Char\naltlocids All alternative location IDs in a DisorderedAtom Array{Char,1}\nx x coordinate of an atom Float64\ny y coordinate of an atom Float64\nz z coordinate of an atom Float64\ncoords coordinates of an atom Array{Float64,1}\noccupancy Occupancy of an atom (default is 1.0) Float64\ntempfactor Temperature factor of an atom (default is 0.0) Float64\nelement Element of an atom (default is \"  \") String\ncharge Charge of an atom (default is \"  \") String\nresidue Residue an atom belongs to Residue\nishetero true if the residue or atom is a hetero residue/atom Bool\nisdisorderedatom true if the atom is disordered Bool\npdbline PDB ATOM/HETATM record for an atom String\nresname Residue name of a residue or atom String\nresnames All residue names in a DisorderedResidue Array{String,1}\nresnumber Residue number of a residue or atom Int\nsequentialresidues Determine if the second residue follows the first in sequence Bool\ninscode Insertion code of a residue or atom Char\nresid Residue ID of an atom or residue (full=true includes chain) String\natomnames Atom names of the atoms in a residue, sorted by serial Array{String,1}\natoms Dictionary of atoms in a residue Dict{String, AbstractAtom}\nisdisorderedres true if the residue has multiple residue names Bool\ndisorderedres Access a particular residue name in a DisorderedResidue Residue\nchain Chain a residue or atom belongs to Chain\nchainid Chain ID of a chain, residue or atom String\nresids Sorted residue IDs in a chain Array{String,1}\nresidues Dictionary of residues in a chain Dict{String, AbstractResidue}\nmodel Model a chain, residue or atom belongs to Model\nmodelnumber Model number of a model, chain, residue or atom Int\nchainids Sorted chain IDs in a model or structure Array{String,1}\nchains Dictionary of chains in a model or structure Dict{String, Chain}\nstructure Structure a model, chain, residue or atom belongs to ProteinStructure\nstructurename Name of the structure an element belongs to String\nmodelnumbers Sorted model numbers in a structure Array{Int,1}\nmodels Dictionary of models in a structure Dict{Int, Model}The strip keyword argument determines whether surrounding whitespace is stripped for atomname, element, charge, resname and atomnames (default true).The coordinates of an atom can be set using x!, y!, z! and coords!."
},

{
    "location": "documentation.html#Manipulating-structures-1",
    "page": "Documentation",
    "title": "Manipulating structures",
    "category": "section",
    "text": "Elements can be looped over to reveal the sub-elements in the correct order:for mod in struc\n    for ch in mod\n        for res in ch\n            for at in res\n                # Do something\n            end\n        end\n    end\nendModels are ordered numerically; chains are ordered by chain ID character ordering, except the empty chain is last; residues are ordered by residue number and insertion code with hetero residues after standard residues; atoms are ordered by atom serial. If you want the first sub-element you can use first. For example first(struc[1]) gets the first chain in model 1. Since the ordering of elements is defined you can use the sort function. For example sort(res) sorts a list of residues as described above, or sort(res, by=resname) will sort them alphabetically by residue name.collect can be used to get arrays of sub-elements. collectatoms, collectresidues, collectchains and collectmodels return arrays of a particular type from a structural element or element array.Selectors are functions passed as additional arguments to these functions. Only elements that return true when passed to all the selector are retained. For example:Command Action Return type\ncollect(struc[\'A\'][50]) Collect the sub-elements of an element, e.g. atoms from a residue Array{AbstractAtom,1}\ncollectresidues(struc) Collect the residues of an element Array{AbstractResidue,1}\ncollectatoms(struc) Collect the atoms of an element Array{AbstractAtom,1}\ncollectatoms(struc, calphaselector) Collect the C-alpha atoms of an element Array{AbstractAtom,1}\ncollectatoms(struc, calphaselector, disorderselector) Collect the disordered C-alpha atoms of an element Array{AbstractAtom,1}The selectors available are:Function Acts on Selects for\nstandardselector AbstractAtom or AbstractResidue Atoms/residues arising from standard (ATOM) records\nheteroselector AbstractAtom or AbstractResidue Atoms/residues arising from hetero (HETATM) records\natomnameselector AbstractAtom Atoms with atom name in a given list\ncalphaselector AbstractAtom C-alpha atoms\ncbetaselector AbstractAtom C-beta atoms, or C-alpha atoms for glycine residues\nbackboneselector AbstractAtom Atoms in the protein backbone (CA, N and C)\nheavyatomselector AbstractAtom Non-hydrogen atoms\nhydrogenselector AbstractAtom Hydrogen atoms\nresnameselector AbstractAtom or AbstractResidue Atoms/residues with residue name in a given list\nwaterselector AbstractAtom or AbstractResidue Atoms/residues with residue name HOH\nnotwaterselector AbstractAtom or AbstractResidue Atoms/residues with residue name not HOH\ndisorderselector AbstractAtom or AbstractResidue Atoms/residues with alternative locationsTo create a new atomnameselector or resnameselector:cdeltaselector(at::AbstractAtom) = atomnameselector(at, [\"CD\"])It is easy to define your own atom, residue, chain or model selectors. The below will collect all atoms with x coordinate less than 0:xselector(at) = x(at) < 0\ncollectatoms(struc, xselector)Alternatively, you can use an anonymous function:collectatoms(struc, at -> x(at) < 0)countatoms, countresidues, countchains and countmodels can be used to count elements with the same selector API. For example:julia> countatoms(struc)\n754\n\njulia> countatoms(struc, calphaselector)\n85\n\njulia> countresidues(struc, standardselector)\n85The sequence of a protein can be retrieved by passing a Chain or array of residues to AminoAcidSequence:julia> AminoAcidSequence(struc[\'A\'], standardselector)\n85aa Amino Acid Sequence:\nRCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRCSee BioSequences.jl for more on how to deal with sequences."
},

{
    "location": "documentation.html#Spatial-calculations-1",
    "page": "Documentation",
    "title": "Spatial calculations",
    "category": "section",
    "text": "Various functions are provided to calculate spatial quantities for proteins:Command Returns\ndistance Minimum distance between two elements\nsqdistance Minimum square distance between two elements\ncoordarray Atomic coordinates in Å of an element as a 2D Array with each column corresponding to one atom\nbondangle Angle between three atoms\ndihedralangle Dihedral angle defined by four atoms\nomegaangle Omega dihedral angle between a residue and the previous residue\nphiangle Phi dihedral angle between a residue and the previous residue\npsiangle Psi dihedral angle between a residue and the next residue\nomegaangles Vector of omega dihedral angles of an element\nphiangles Vector of phi dihedral angles of an element\npsiangles Vector of psi dihedral angles of an element\nramachandranangles Vectors of phi and psi angles of an element\nContactMap ContactMap of two elements, or one element with itself\nDistanceMap DistanceMap of two elements, or one element with itself\nshowcontactmap Print a representation of a ContactMap to stdout or a specified IO instance\nrmsd RMSD between two elements of the same size - assumes they are superimposed\ndisplacements Vector of displacements between two elements of the same size - assumes they are superimposedThe omegaangle, phiangle and psiangle functions can take either a pair of residues or a chain and a position. The omegaangle and phiangle functions measure the angle between the residue at the given index and the one before. The psiangle function measures between the given index and the one after.For example:julia> distance(struc[\'A\'][10], struc[\'A\'][20])\n10.782158874733762\n\njulia> rad2deg(bondangle(struc[\'A\'][50][\"N\"], struc[\'A\'][50][\"CA\"], struc[\'A\'][50][\"C\"]))\n110.77765846083398\n\njulia> rad2deg(dihedralangle(struc[\'A\'][50][\"N\"], struc[\'A\'][50][\"CA\"], struc[\'A\'][50][\"C\"], struc[\'A\'][51][\"N\"]))\n-177.38288114072924\n\njulia> rad2deg(psiangle(struc[\'A\'][50], struc[\'A\'][51]))\n-177.38288114072924\n\njulia> rad2deg(psiangle(struc[\'A\'], 50))\n-177.38288114072924ContactMap takes in a structural element or a list, such as a Chain or Vector{Atom}, and returns a ContactMap object showing the contacts between the elements for a specified distance. ContactMap can also be given two structural elements as arguments, in which case a non-symmetrical 2D array is returned showing contacts between the elements. The underlying BitArray{2} for ContactMap contacts can be accessed with contacts.data if required.julia> contacts = ContactMap(collectatoms(struc[\'A\'], cbetaselector), 8.0)\nContact map of size (85, 85)A plot recipe is defined for this so it can shown with Plots.jl:using Plots\nplot(contacts)(Image: contactmap)For a quick, text-based representation of a ContactMap use showcontactmap.DistanceMap works in an analogous way to ContactMap and gives a map of the distances. It can also be plotted:dists = DistanceMap(collectatoms(struc[\'A\'], cbetaselector))\nusing Plots\nplot(dists)(Image: distancemap)"
},

{
    "location": "documentation.html#Downloading-PDB-files-1",
    "page": "Documentation",
    "title": "Downloading PDB files",
    "category": "section",
    "text": "To download a PDB file to a specified directory:downloadpdb(\"1EN2\", pdb_dir=\"path/to/pdb/directory\")To download multiple PDB files to a specified directory:downloadpdb([\"1EN2\", \"1ALW\", \"1AKE\"], pdb_dir=\"path/to/pdb/directory\")To download a PDB file in PDB, XML, MMCIF or MMTF format use the file_format argument:downloadpdb(\"1ALW\", pdb_dir=\"path/to/pdb/directory\", file_format=MMTF)\n\n# To get XML\ndownloadpdb(\"1ALW\", pdb_dir=\"path/to/pdb/directory\", file_format=PDBXML)To apply a function to a downloaded file and delete the file afterwards:downloadpdb(f, \"1ALW\")Or, using Julia\'s do syntax:downloadpdb(\"1ALW\") do fp\n    s = read(fp, PDB)\n    # Do something\nend"
},

{
    "location": "documentation.html#Reading-PDB-files-1",
    "page": "Documentation",
    "title": "Reading PDB files",
    "category": "section",
    "text": "To parse an existing PDB file into a Structure-Model-Chain-Residue-Atom framework:julia> struc = read(\"/path/to/pdb/file.pdb\", PDB)\nProteinStructure 1EN2.pdb with 1 models, 1 chains (A), 85 residues, 754 atomsRead a mmCIF file instead by replacing PDB with MMCIF. Various options can be set through optional keyword arguments when parsing PDB/mmCIF files:Keyword Argument Description\nstructure_name::AbstractString The name given to the returned ProteinStructure; defaults to the file name\nremove_disorder::Bool=false Whether to remove atoms with alt loc ID not \' \' or \'A\'.\nread_std_atoms::Bool=true Whether to read standard ATOM records.\nread_het_atoms::Bool=true Whether to read HETATOM records.The function readpdb provides an alternative way to read PDB files with a similar interface to downloadpdb. To parse a PDB file by specifying the PDB ID and PDB directory:struc = readpdb(\"1EN2\", pdb_dir=\"/path/to/pdb/directory\")The same keyword arguments are taken as read above, plus pdb_dir and ba_number.Use retrievepdb to download and parse a PDB file into a Structure-Model-Chain-Residue-Atom framework in a single line:julia> struc = retrievepdb(\"1ALW\", pdb_dir=\"path/to/pdb/directory\")\nINFO: Downloading PDB: 1ALW\nProteinStructure 1ALW.pdb with 1 models, 2 chains (A,B), 346 residues, 2928 atoms"
},

{
    "location": "documentation.html#Writing-PDB-files-1",
    "page": "Documentation",
    "title": "Writing PDB files",
    "category": "section",
    "text": "PDB format files can be written:writepdb(\"1EN2_out.pdb\", struc)Any element type can be given as input to writepdb. Atom selectors can also be given as additional arguments:# Only backbone atoms are written out\nwritepdb(\"1EN2_out.pdb\", struc, backboneselector)The first argument can also be a stream. To write mmCIF format files, use the writemmcif function with similar arguments. A MMCIFDict can also be written using writemmcif:writemmcif(\"1EN2_out.dic\", mmcif_dict)Multi-character chain IDs can be written to mmCIF files but will throw an error when written to a PDB file as the PDB file format only has one character allocated to the chain ID.If you want the PDB record line for an Atom, use pdbline. For example:julia> pdbline(at)\n\"HETATM  101  C  A  X B  20      10.500  20.123  -5.123  0.50 50.13           C1+\"If you want to generate a PDB record line from values directly, do so using an AtomRecord:julia> pdbline(AtomRecord(false, 669, \"CA\", \' \', \"ILE\", \"A\", 90, \' \', [31.743, 33.11, 31.221], 1.00, 25.76, \"C\", \"\"))\n\"ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  \"This can be useful when writing PDB files from your own data structures."
},

{
    "location": "documentation.html#RCSB-PDB-utility-functions-1",
    "page": "Documentation",
    "title": "RCSB PDB utility functions",
    "category": "section",
    "text": "To get the list of all PDB entries:l = pdbentrylist()To download the entire RCSB PDB database in your preferred file format:downloadentirepdb(pdb_dir=\"path/to/pdb/directory\", file_format=MMTF)This operation takes a lot of disk space and time to complete (depending on internet connection).To update your local PDB directory based on the weekly status list of new, modified and obsolete PDB files from the RCSB server:updatelocalpdb(pdb_dir=\"path/to/pdb/directory\", file_format=MMTF)Obsolete PDB files are stored in the auto-generated obsolete directory inside the specified local PDB directory.To maintain a local copy of the entire RCSB PDB database, run the downloadentirepdb function once to download all PDB files and set up a CRON job or similar to run updatelocalpdb function once a week to keep the local PDB directory up to date with the RCSB server.There are a few more functions that may be useful:Function Returns Return type\npdbentrylist List of all PDB entries from the RCSB server Array{String,1}\npdbstatuslist List of PDB entries from a specified RCSB weekly status list URL Array{String,1}\npdbrecentchanges Added, modified and obsolete PDB lists from the recent RCSB weekly status files Tuple{Array{String,1},Array{String,1}, Array{String,1}}\npdbobsoletelist List of all obsolete PDB entries Array{String,1}\ndownloadallobsoletepdb Downloads all obsolete PDB files from the RCSB PDB server Array{String,1}"
},

{
    "location": "documentation.html#Visualising-structures-1",
    "page": "Documentation",
    "title": "Visualising structures",
    "category": "section",
    "text": "The Bio3DView.jl package can be used to visualise molecular structures. For example:using Bio3DView\nusing Blink\nviewpdb(\"1CRN\")(Image: viewpdb)struc = retrievepdb(\"1AKE\")\nviewstruc(struc[\'A\'], surface=Surface(Dict(\"colorscheme\"=> \"greenCarbon\")))(Image: viewstruc)Here they are shown as static images but they are interactive when using Bio3DView.jl. See the Bio3DView.jl tutorial for more information."
},

{
    "location": "documentation.html#Examples-1",
    "page": "Documentation",
    "title": "Examples",
    "category": "section",
    "text": "A few further examples of BioStructures usage are given below.A) Plot the temperature factors of a protein:using Plots\ncalphas = collectatoms(struc, calphaselector)\nplot(resnumber.(calphas),\n     tempfactor.(calphas),\n     xlabel=\"Residue number\",\n     ylabel=\"Temperature factor\",\n     label=\"\")B) Print the PDB records for all C-alpha atoms within 5 Angstrom of residue 38:for at in calphas\n    if distance(struc[\'A\'][38], at) < 5.0 && resnumber(at) != 38\n        println(pdbline(at))\n    end\nendC) Find the residues at the interface of a protein-protein interaction:for res_a in collectresidues(struc[\"A\"], standardselector)\n    for res_b in collectresidues(struc[\"B\"], standardselector)\n        if distance(res_a, res_b) < 5.0\n            println(resnumber(res_a), \"A \", resnumber(res_b), \"B\")\n        end\n    end\nendD) Show the Ramachandran phi/psi angle plot of a structure:using Plots\nphi_angles, psi_angles = ramachandranangles(struc, standardselector)\nscatter(rad2deg.(phi_angles),\n     rad2deg.(psi_angles),\n     title=\"Ramachandran plot\",\n     xlabel=\"Phi / degrees\",\n     ylabel=\"Psi / degrees\",\n     label=\"\",\n     xticks=[-180, -90, 0, 90, 180],\n     yticks=[-180, -90, 0, 90, 180],\n     xlims=(-180, 180),\n     ylims=(-180, 180))E) Calculate the RMSD and displacements between the heavy (non-hydrogen) atoms of two models in an NMR structure:downloadpdb(\"1SSU\")\nstruc_nmr = read(\"1SSU.pdb\", PDB)\nrmsd(struc_nmr[5], struc_nmr[10], heavyatomselector)\ndisplacements(struc_nmr[5], struc_nmr[10], heavyatomselector)F) Calculate the cysteine fraction of every structure in the PDB:l = pdbentrylist()\nfor p in l\n    downloadpdb(p, file_format=MMCIF) do fp\n        s = read(fp, MMCIF)\n        nres = countresidues(s, standardselector)\n        if nres > 0\n            frac = countresidues(s, standardselector, x -> resname(x) == \"CYS\") / nres\n            println(p, \"  \", round(frac, digits=2))\n        end\n    end\nendG) Interoperability is possible with other packages in the Julia ecosystem. For example, use NearestNeighbors.jl to find the 10 nearest residues to each residue:using NearestNeighbors\nstruc = retrievepdb(\"1AKE\")\nca = coordarray(struc[\"A\"], cbetaselector)\nkdtree = KDTree(ca; leafsize=10)\nidxs, dists = knn(kdtree, ca, 10, true)"
},

]}
