var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = CamiXon","category":"page"},{"location":"#CamiXon.jl","page":"Home","title":"CamiXon.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package for image analysis of backscattered light","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Finite-difference-methods","page":"Home","title":"Finite-difference methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"f_diff_weight(k::Int, i::Int)\nf_diff_weights(k::Int)\nf_diff_weights_array(kmax::Int)\nf_diff_expansion_weights(coeffs, ∇)\nf_diff_expansion_coeffs_interpolation(k::Int, x::T) where T<:Real","category":"page"},{"location":"#CamiXon.f_diff_weight-Tuple{Int64, Int64}","page":"Home","title":"CamiXon.f_diff_weight","text":"f_diff_weight(k::Int, i::Int)\n\nWeight coefficient\n\nc_j^k=(-1)^jbinomkj\n\nof the k^th-order finite difference operator nabla^k and corresponding to the function value fn-j.\n\nExample:\n\nk = 5; i = 3\nf_diff_weight(k, i)\n -10\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.f_diff_weights-Tuple{Int64}","page":"Home","title":"CamiXon.f_diff_weights","text":"f_diff_weights(k::Int)\n\nWeight coefficients c_0^k ldots c_k^k defining the k^th-order finite difference operator,\n\nnabla^k fn = sum_j=0^k c_i^kfn-j\n\nwhere fn fn-k are elements of a tabulated analytic function.\n\nExample:\n\nk = 3\nf_diff_weights(k)\n4-element Vector{Int64}:\n  1\n -3\n  3\n -1\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.f_diff_weights_array-Tuple{Int64}","page":"Home","title":"CamiXon.f_diff_weights_array","text":"f_diff_weights_array(kmax::Int)\n\nCollection of weight coefficients c_0^k ldots c_j^k defining the finite difference operator nabla^j  (0le jle k).\n\nExample:\n\nkmax = 3\n∇ = f_diff_weights_array(kmax)\n4-element Vector{Vector{Int64}}:\n [1]\n [1, -1]\n [1, -2, 1]\n [1, -3, 3, -1]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.f_diff_expansion_weights-Tuple{Any, Any}","page":"Home","title":"CamiXon.f_diff_expansion_weights","text":"f_diff_expansion_weights(a, ∇)\n\nSummation weights b_0^k ldots b_k^k corresponding to the expansion coefficients a_0^k ldots a_k^k of a k^th-order finite-difference expansion,\n\nsum_p=0^ka_pnabla^pfn=sum_j=0^kb_j^kfn-j\n\nwhere fn fn-k are elements of a tabulated analytic function.\n\nExample:\n\nk=5\n∇ = f_diff_weights_array(k)\na = UnitRange(0,k)\nb = f_diff_expansion_weights(a, ∇)\n6-element Vector{Int64}:\n  15\n -55\n  85\n -69\n  29\n  -5\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.f_diff_expansion_coeffs_interpolation-Union{Tuple{T}, Tuple{Int64, T}} where T<:Real","page":"Home","title":"CamiXon.f_diff_expansion_coeffs_interpolation","text":"f_diff_expansion_coeffs_interpolation(k::Int, x::T) where T<:Real\n\nFinite-difference expansion coefficients l_p(x) for lagrangian interpolation of a tabulated analytic function at offset position 0le xle -k,\n\nfn+x =sum_p=0^kl_p(x)nabla^pfn = sum_j=0^kr_j^k(x)fn-j\n\nwhere l_0equiv 0 and\n\nl_p(x) = x(x+1)(x+2)cdots(x+p-1)p\n\nwith p=1 ldots k, and fn fn-k are elements of the function table. The lagrangian interpolation weights r_0^k ldots r_k^k are calculated with the function r = f_diff_expansion_weights(l, ∇).\n\n\n\nk=3\n∇ = f_diff_weights_array(k)\nx=-1\nl = f_diff_expansion_coeffs_interpolation(k,x)\nr = f_diff_expansion_weights(l, ∇)\nprintln(l,r)\n [1, -1, 0, 0][0, 1, 0, 0]\n\n\n\n\n\n","category":"method"},{"location":"#FITS","page":"Home","title":"FITS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FITS stands for 'Flexible Image Transport System'. This is an open standard origionally developed for the astronomy community to store telescope images together with tables of spectral information. Over the years it has developed into a scientific standard - http://fits.gsfc.nasa.gov/iaufwg.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Within CamiXion only the basic FITS functionality is implemented for users not requiring celestal coordinates. The user can create, read and extend .fits files as well as create, edit and delete user-defined metainformation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A FITS file consists of a sequence of one or more header-data-units (HDUs), each containing a data block preceeded by header records of metainformation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"By the command f = fits_read(filnam) we asign a collection of FITS_HDU objects from the file filnam to the variable f.","category":"page"},{"location":"#FITS-Types","page":"Home","title":"FITS - Types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FITS_HDU\nFITS_header\nFITS_data\nFITS_table\nFITS_name","category":"page"},{"location":"#CamiXon.FITS_HDU","page":"Home","title":"CamiXon.FITS_HDU","text":"FITS_HDU\n\nObject to hold a single \"Header-Data Unit\" (HDU).\n\nThe fields are\n\n.filename::String: name of the corresponding FITS file\n.hduindex::Int: identifier (a file may contain more than one HDU)\n.header::T, where T=FITS_header: the header object\n.dataobject::V, where V=FITS_data: the data object\n\n\n\n\n\n","category":"type"},{"location":"#CamiXon.FITS_header","page":"Home","title":"CamiXon.FITS_header","text":"FITS_header\n\nObject to hold the header information of a FITS_HDU.\n\nThe fields are:\n\n.hduindex::Int: identifier (a file may contain more than one HDU)\n.records::Array{String,1}:  the header formated as an array of strings of 80 ASCII characters\n.keys::Array{String,1}: keys[i] - key corresponding to records[i] (record of index i)\n.values::Array{Any,1}: value[i] - corresponding to records[i]\n.comments: comments[i] - comment corresponding to records[i]\n.dict::Dict{String,Any}: Dictionary key[i] => value[i]\n.maps::Dict{String,Int}: Dictionary key[i] => i\n\n\n\n\n\n","category":"type"},{"location":"#CamiXon.FITS_data","page":"Home","title":"CamiXon.FITS_data","text":"FITS_data\n\nObject to hold the data of the FITS_HDU of given hduindex and hdutype.\n\nThe fields are:\n\n.hduindex::Int: identifier (a file may contain more than one HDU)\n.hdutype::String: accepted types are 'PRIMARY', 'IMAGE' and 'TABLE'\n.data::Any: in the from appropriate for the hdutype\n\n\n\n\n\n","category":"type"},{"location":"#CamiXon.FITS_table","page":"Home","title":"CamiXon.FITS_table","text":"FITS_table\n\nObject to hold the data of a TABLE HDU (a FITS_HDU for ASCII tables). It contains the data in the form of records (rows) of ASCII strings.\n\nThe fields are:\n\n.hduindex::Int: identifier (a file may contain more than one HDU)\n.rows::Array{String,1}: the table formated as an array of rows of ASCII strings\n\n\n\n\n\n","category":"type"},{"location":"#CamiXon.FITS_name","page":"Home","title":"CamiXon.FITS_name","text":"FITS_name\n\nFITS object to decompose the names of .fits files.\n\nThe fields are:\n\n.name::String: for filename 'p#.fits' this is 'p#.fits'\n.prefix::String: for filename 'p#.fits' this is 'p'\n.numerator::String: for filename 'p#.fits' this is '#', a serial number (e.g., '3') or a range (e.g., '3-7')\n.extension::String:  for filename 'p#.fits' this is '.fits'.\n\n\n\n\n\n","category":"type"},{"location":"#FITS-HDU-Methods","page":"Home","title":"FITS - HDU Methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"fits_info(hdu::FITS_HDU)\nparse_FITS_TABLE(hdu::FITS_HDU)","category":"page"},{"location":"#CamiXon.fits_info-Tuple{FITS_HDU}","page":"Home","title":"CamiXon.fits_info","text":"fits_info(hdu)\n\nPrint metafinformation and data of given FITS_HDU\n\nExample:\n\nstrExample = \"remove.fits\"\ndata = [11,21,31,12,22,23,13,23,33]\ndata = reshape(data,(3,3,1))\nfits_create(strExample, data; protect=false)\n\nf = fits_read(strExample)\nfits_info(f[1])\n\n  File: remove.fits\n  hdu: 1\n  hdutype: PRIMARY\n  DataType: Int64\n  Datasize: (3, 3, 1)\n\n  Metainformation:\n  SIMPLE  =                    T / file does conform to FITS standard\n  BITPIX  =                   64 / number of bits per data pixel\n  NAXIS   =                    3 / number of data axes\n  NAXIS1  =                    3 / length of data axis 1\n  NAXIS2  =                    3 / length of data axis 2\n  NAXIS3  =                    1 / length of data axis 3\n  BZERO   =                  0.0 / offset data range to that of unsigned integer\n  BSCALE  =                  1.0 / default scaling factor\n  EXTEND  =                    T / FITS dataset may contain extensions\n  COMMENT    Primary FITS HDU    / http://fits.gsfc.nasa.gov/iaufwg\n  END\n\n  3×3×1 Array{Int64, 3}:\n  [:, :, 1] =\n   11  12  13\n   21  22  23\n   31  23  33\n\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.parse_FITS_TABLE-Tuple{FITS_HDU}","page":"Home","title":"CamiXon.parse_FITS_TABLE","text":"parse_FITS_TABLE(hdu)\n\nParse FITS_TABLE (ASCII table) into a Vector of its columns for further processing by the user. Default formatting in ISO 2004 FORTRAN data format specified by keys \"TFORMS1\" - \"TFORMSn\") Display formatting in ISO 2004 FORTRAN data format (\"TDISP1\" - \"TDISPn\") prepared for user editing.\n\nExample:\n\nstrExample = \"example.fits\"\ndata = [10, 20, 30]\nfits_create(strExample, data; protect=false)\n\nt1 = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]\nt2 = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]\nt3 = [1.23,2.12,3.,4.,5.]\nt4 = ['a','b','c','d','e']\nt5 = [\"a\",\"bb\",\"ccc\",\"dddd\",\"ABCeeaeeEEEEEEEEEEEE\"]\ndata = [t1,t2,t3,t4,t5]\nfits_extend(strExample, data, \"TABLE\")\n\nf = fits_read(strExample)\nd = f[2].header.dict\nd = [get(d,\"TFORM$i\",0) for i=1:5]; println(strip.(d))\n  SubString{String}[\"'E6.1    '\", \"'I4      '\", \"'F4.2    '\", \"'A1      '\", \"'A20     '\"]\n\nf[2].dataobject.data                            # this is the table hdu\n  5-element Vector{String}:\n   \"1.0e-6 1086 1.23 a a                    \"\n   \"2.0e-6 1036 2.12 b bb                   \"\n   \"3.0e-6 1055 3.0  c ccc                  \"\n   \"4.0e-6 1070 4.0  d dddd                 \"\n   \"5.0e-6 1071 5.0  e ABCeeaeeEEEEEEEEEEEE \"\n\nparse_FITS_TABLE(f[2])\n  5-element Vector{Vector{T} where T}:\n   [1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6, 5.0e-6]\n   [1086, 1036, 1055, 1070, 1071]\n   [1.23, 2.12, 3.0, 4.0, 5.0]\n   [\"a\", \"b\", \"c\", \"d\", \"e\"]\n   [\"a                   \", \"bb                  \", \"ccc                 \", \"dddd                \", \"ABCeeaeeEEEEEEEEEEEE\"]\n\n\n\n\n\n","category":"method"},{"location":"#FITS-File-Methods","page":"Home","title":"FITS - File Methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"cast_FITS_name(filename::String)\nfits_combine(filnamFirst::String, filnamLast::String; protect=true)\nfits_copy(filenameA::String, filenameB::String=\" \"; protect=true)\nfits_create(filename::String, data=[]; protect=true)\nfits_extend(filename::String, data_extend, hdutype=\"IMAGE\")\nfits_read(filename::String)","category":"page"},{"location":"#CamiXon.cast_FITS_name-Tuple{String}","page":"Home","title":"CamiXon.cast_FITS_name","text":"cast_FITS_name(filename::String)\n\nDecompose the FITS filename 'filnam.fits' into its name, prefix, numerator and extension.\n\nExamples:\n\nstrExample = \"T23.01.fits\"\nf = cast_FITS_name(strExample)\nFITS_name(\"T23.01\", \"T23.\", \"01\", \".fits\")\n\nf.name, f.prefix, f.numerator, f.extension\n(\"T23.01\", \"T23.\", \"01\", \".fits\")\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.fits_combine-Tuple{String, String}","page":"Home","title":"CamiXon.fits_combine","text":"fits_combine(strFirst, strLast [; protect=true])\n\nCopy \"filenameFirst\" to \"filenameLast\" (with mandatory \".fits\" extension)\n\nKey:\n\nprotect::Bool: overwrite protection\n\nExample:\n\nfits_combine(\"T01.fits\", \"T22.fits\")\n  'T01-T22.fits': file created\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.fits_copy","page":"Home","title":"CamiXon.fits_copy","text":"fits_copy(filenameA [, filenameB=\"\" [; protect=true]])\n\nCopy \"filenameA\" to \"filenameB\" (with mandatory \".fits\" extension) Key:\n\nprotect::Bool: overwrite protection\n\nExamples:\n\nfits_copy(\"T01.fits\")\n  'T01.fits' was saved as 'T01 - Copy.fits'\n\nfits_copy(\"T01.fits\", \"T01a.fits\")\n  FitsError: 'T01a.fits' in use (set ';protect=false' to lift overwrite protection)\n\nfits_copy(\"T01.fits\", \"T01a.fits\"; protect=false)\n  'T01.fits' was saved as 'T01a.fits'\n\n\n\n\n\n","category":"function"},{"location":"#CamiXon.fits_create","page":"Home","title":"CamiXon.fits_create","text":"fits_create(filename [, data [; protect=true]])\n\nCreate FITS file of given filename [, optional data block [, default overwrite protection]] and return Array of HDUs. Key:\n\nprotect::Bool: overwrite protection\n\nExamples:\n\nstrExample = \"minimal.fits\"\nfits_create(strExample;protect=false)\n\nf = fits_read(strExample)\na = f[1].dataobject.data\nb = f[1].header.keys\nprintln(a);println(b)\n  Any[]\n  [\"SIMPLE\", \"NAXIS\", \"EXTEND\", \"COMMENT\", \"END\"]\n\nstrExample = \"remove.fits\"\ndata = [11,21,31,12,22,23,13,23,33]\ndata = reshape(data,(3,3,1))\nfits_create(strExample, data; protect=false)\n\nf = fits_read(strExample)\nfits_info(f[1])\n\n  File: remove.fits\n  hdu: 1\n  hdutype: PRIMARY\n  DataType: Int64\n  Datasize: (3, 3, 1)\n\n  Metainformation:\n  SIMPLE  =                    T / file does conform to FITS standard\n  BITPIX  =                   64 / number of bits per data pixel\n  NAXIS   =                    3 / number of data axes\n  NAXIS1  =                    3 / length of data axis 1\n  NAXIS2  =                    3 / length of data axis 2\n  NAXIS3  =                    1 / length of data axis 3\n  BZERO   =                  0.0 / offset data range to that of unsigned integer\n  BSCALE  =                  1.0 / default scaling factor\n  EXTEND  =                    T / FITS dataset may contain extensions\n  COMMENT    Primary FITS HDU    / http://fits.gsfc.nasa.gov/iaufwg\n  END\n\n  3×3×1 Array{Int64, 3}:\n  [:, :, 1] =\n   11  12  13\n   21  22  23\n   31  23  33\n\n\n\n\n\n","category":"function"},{"location":"#CamiXon.fits_extend","page":"Home","title":"CamiXon.fits_extend","text":"fits_extend(filename, data_extend, hdutype=\"IMAGE\")\n\nExtend the FITS file of given filename with the data of hdutype from data_extend  and return Array of HDUs.\n\nExamples:\n\nstrExample = \"test_example.fits\"\ndata = [0x0000043e, 0x0000040c, 0x0000041f]\nfits_create(strExample, data, protect=false)\n\nf = fits_read(strExample)\na = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]\nb = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]\nc = [1.23,2.12,3.,4.,5.]\nd = ['a','b','c','d','e']\ne = [\"a\",\"bb\",\"ccc\",\"dddd\",\"ABCeeaeeEEEEEEEEEEEE\"]\ndata = [a,b,c,d,e]\nfits_extend(strExample, data, \"TABLE\")\n\nf = fits_read(strExample)\nf[2].dataobject.data\n  5-element Vector{String}:\n   \"1.0e-6 1086 1.23 a a                    \"\n   \"2.0e-6 1036 2.12 b bb                   \"\n   \"3.0e-6 1055 3.0  c ccc                  \"\n   \"4.0e-6 1070 4.0  d dddd                 \"\n   \"5.0e-6 1071 5.0  e ABCeeaeeEEEEEEEEEEEE \"\n\nrm(strExample); f = data = a = b = c = d = e = nothing\n\n\n\n\n\n","category":"function"},{"location":"#CamiXon.fits_read-Tuple{String}","page":"Home","title":"CamiXon.fits_read","text":"fits_read(filename)\n\nRead FITS file and return Array of FITS_HDUs\n\nExample:\n\nstrExample = \"minimal.fits\"\nfits_create(strExample;protect=false)\n\nf = fits_read(strExample)\nf[1].dataobject.data\n  Any[]\n\nrm(strExample); f = nothing\n\n\n\n\n\n","category":"method"},{"location":"#FITS-Key-Methods","page":"Home","title":"FITS - Key Methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"fits_add_key(filename::String, hduindex::Int, key::String, val::Real, com::String)\nfits_delete_key(filename::String, hduindex::Int, key::String)\nfits_edit_key(filename::String, hduindex::Int, key::String, val::Real, com::String)\nfits_rename_key(filename::String, hduindex::Int, keyold::String, keynew::String)","category":"page"},{"location":"#CamiXon.fits_add_key-Tuple{String, Int64, String, Real, String}","page":"Home","title":"CamiXon.fits_add_key","text":"fits_add_key(filename, hduindex, key, value, comment)\n\nAdd a header record of given 'key, value and comment' to 'HDU[hduindex]' of file with name 'filename'\n\nExample:\n\nstrExample=\"minimal.fits\"\nfits_create(strExample;protect=false)\nfits_add_key(strExample, 1, \"KEYNEW1\", true, \"FITS dataset may contain extension\")\n\nf = fits_read(strExample)\nfits_info(f[1])\n\n  File: minimal.fits\n  hdu: 1\n  hdutype: PRIMARY\n  DataType: Any\n  Datasize: (0,)\n\n  Metainformation:\n  SIMPLE  =                    T / file does conform to FITS standard\n  NAXIS   =                    0 / number of data axes\n  EXTEND  =                    T / FITS dataset may contain extensions\n  COMMENT    Primary FITS HDU    / http://fits.gsfc.nasa.gov/iaufwg\n  KEYNEW1 =                    T / FITS dataset may contain extension\n  END\n\n  Any[]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.fits_delete_key-Tuple{String, Int64, String}","page":"Home","title":"CamiXon.fits_delete_key","text":"fits_delete_key(filename, hduindex, key)\n\nDelete a header record of given key, value and comment to FITS_HDU[hduindex] of file with name  'filename'\n\nExamples:\n\nstrExample=\"minimal.fits\"\nfits_create(strExample;protect=false)\nfits_add_key(strExample, 1, \"KEYNEW1\", true, \"this is record 5\")\n\nf = fits_read(strExample)\nget(f[1].header.maps,\"KEYNEW1\",0)\n  5\n\nfits_delete_key(strExample, 1, \"KEYNEW1\")\n\nf = fits_read(strExample)\nget(f[1].header.maps,\"KEYNEW1\",0)\n  0\n\nfits_delete_key(filnam, 1, \"NAXIS\")\n 'NAXIS': cannot be deleted (key protected under FITS standard)\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.fits_edit_key-Tuple{String, Int64, String, Real, String}","page":"Home","title":"CamiXon.fits_edit_key","text":"fits_edit_key(filename, hduindex, key, value, comment)\n\nEdit a header record of given 'key, value and comment' to 'HDU[hduindex]' of file with name 'filename'\n\nExample:\n\ndata = DateTime(\"2020-01-01T00:00:00.000\")\nstrExample=\"minimal.fits\"\nfits_create(strExample;protect=false)\nfits_add_key(strExample, 1, \"KEYNEW1\", true, \"this is record 5\")\nfits_edit_key(strExample, 1, \"KEYNEW1\", data, \"record 5 changed to a DateTime type\")\n\nf = fits_read(strExample)\nfits_info(f[1])\n\n  File: minimal.fits\n  hdu: 1\n  hdutype: PRIMARY\n  DataType: Any\n  Datasize: (0,)\n\n  Metainformation:\n  SIMPLE  =                    T / file does conform to FITS standard\n  NAXIS   =                    0 / number of data axes\n  EXTEND  =                    T / FITS dataset may contain extensions\n  COMMENT    Primary FITS HDU    / http://fits.gsfc.nasa.gov/iaufwg\n  KEYNEW1 = '2020-01-01T00:00:00' / record 5 changed to a DateTime type\n  END\n\n  Any[]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.fits_rename_key-Tuple{String, Int64, String, String}","page":"Home","title":"CamiXon.fits_rename_key","text":"fits_rename_key(filename, hduindex, keyold, kewnew)\n\nRename the key of a header record of file with name 'filename'\n\nExample:\n\nstrExample=\"minimal.fits\"\nfits_create(strExample;protect=false)\nfits_add_key(strExample, 1, \"KEYNEW1\", true, \"this is record 5\")\nfits_rename_key(strExample, 1, \"KEYNEW1\",  \"KEYNEW2\")\n\nf = fits_read(strExample)\nfits_info(f[1])\n\n  File: minimal.fits\n  hdu: 1\n  hdutype: PRIMARY\n  DataType: Any\n  Datasize: (0,)\n\n  Metainformation:\n  SIMPLE  =                    T / file does conform to FITS standard\n  NAXIS   =                    0 / number of data axes\n  EXTEND  =                    T / FITS dataset may contain extensions\n  COMMENT    Primary FITS HDU    / http://fits.gsfc.nasa.gov/iaufwg\n  KEYNEW2 =                    T / this is record 5\n  END\n\n  Any[]\n\n\n\n\n\n","category":"method"},{"location":"#FORTRAN","page":"Home","title":"FORTRAN","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FORTRAN_format\ncast_FORTRAN_format(str::String)\ncast_FORTRAN_datatype(str::String)","category":"page"},{"location":"#CamiXon.FORTRAN_format","page":"Home","title":"CamiXon.FORTRAN_format","text":"FORTRAN_format\n\nObject to hold a FORTRAN format specifier decomposed into its fields. Accepted datatype specifiers are:  Aw,  Iw,  Fw.d,  Ew.d,  Dw.d. Accepted output formating specifiers are: Aw,  Iw.m,  Bw.m,  Ow.m,  Zw.m,  Fw.d,  Ew.dEe,  ENw.d,  ESw.d,  Gw.dEe,  Dw.dEe. Notation: 'w' - width, 'm' (optional) - minimum number of digits, 'd' - number of digits to right of decimal, 'e' - number of digits in exponent; 'N'/'S' (optional) - switch for engineering/scientific formating of the 'E' type.\n\nThe fields are:\n\nType::String: primary FORTRAN datatype\nTypeChar::Char: primary FORTRAN datatype character\nEngSci::Union{Char,Nothing}: secundary datatype character (N for engineering, S for scientific)\nwidth::Int: width of numeric field\nnmin::Int: minimum number of digits displayed\nndec::Int: number of digits to right of decimal\nnexp::Int: number of digits in exponent\n\n\n\n\n\n","category":"type"},{"location":"#CamiXon.cast_FORTRAN_format-Tuple{String}","page":"Home","title":"CamiXon.cast_FORTRAN_format","text":"cast_FORTRAN_format(format::String)\n\nDecompose the format specifier 'format' into its fields and cast this into the FORTRAN_format object. Allowed format specifiers are of the types: Aw, Iw.m, Bw.m, Ow.m, Zw.m, Fw.d, Ew.dEe, ENw.d, ESw.d, Gw.dEe, Dw.dEe, with: 'w' - width, 'm' (optional) - minimum number of digits, 'd' - number of digits to right of decimal, 'e' - number of digits in exponent; 'N'/'S' (optional) - switch for engineering/scientific formating of the 'E' type.\n\nExamples:\n\nt = cast_FORTRAN_format(\"I10\")\nFORTRAN_format(\"Iw\", 'I', nothing, 10, 0, 0, 0)\n\nt = cast_FORTRAN_format(\"I10.12\")\nFORTRAN_format(\"Iw.m\", 'I', nothing, 10, 12, 0, 0)\n\nt = cast_FORTRAN_format(\"E10.5E3\")\nFORTRAN_format(\"Ew.dEe\", 'E', nothing, 10, 0, 5, 3)\n\nt.Type, t.TypeChar, t.EngSci, t.width, t.nmin, t.ndec, t.nexp\n(\"Ew.dEe\", 'E', nothing, 10, 0, 5, 3)\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.cast_FORTRAN_datatype-Tuple{String}","page":"Home","title":"CamiXon.cast_FORTRAN_datatype","text":"cast_FORTRAN_datatype(format::String)\n\nDecompose the format specifier 'format' into its fields and cast this into the FORTRAN_format object. Allowed format specifiers are of the types: Aw, Iw, Fw.d, Ew.d, Dw.d, where: 'w' - width, 'd' - number of digits to right of decimal point.\n\nExamples:\n\nf = cast_FORTRAN_datatype(\"I10\")\nFORTRAN_format(\"Iw\", 'I', nothing, 10, 0, 0, 0)\n\nf = cast_FORTRAN_datatypet(\"F10.4\")\nFORTRAN_format(\"Fw.d\", 'F', nothing, 10, 0, 4, 0)\n\nf = cast_FORTRAN_datatype(\"E10.5\")\nFORTRAN_format(\"Ew.d\", 'E', nothing, 10, 0, 5, 0)\n\nt.Type, t.TypeChar, t.EngSci, t.width, t.nmin, t.ndec, t.nexp\n(\"Ew.d\", 'E', nothing, 10, 0, 5, 0)\n\n\n\n\n\n","category":"method"},{"location":"#Plotting","page":"Home","title":"Plotting","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"step125(x::Real)\nselect125(x)\nsteps(x::Vector{T} where T<:Real)\nstepcenters(x::Vector{T} where T<:Real)\nstepedges(x::Vector{T} where T<:Real)\nedges(px, Δx=1.0, x0=0.0)","category":"page"},{"location":"#CamiXon.step125-Tuple{Real}","page":"Home","title":"CamiXon.step125","text":"step125(x)\n\nStep used for deviding the number x in steps according to 1-2-5 scheme\n\nExamples:\n\nstep125.([5,10,21.3,50,100.1])\n5-element Vector{Int64}:\n  1\n  2\n  5\n 10\n 20\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.select125-Tuple{Any}","page":"Home","title":"CamiXon.select125","text":"select125(x)\n\nSelect elements of the collection x by index according to 1-2-5 scheme\n\nExamples:\n\nx = [1,2,4,6,8,10,13,16,18,20,40,60,80,100]\nselect125(x)\n [2, 6, 10, 16, 20, 60, 100]\n\nx = string.(x)\nselect125(x)\n [\"2\", \"6\", \"10\", \"16\", \"20\", \"60\", \"100\"]\n\nx = 1:100\nselect125(x)\n [20, 40, 60, 80, 100]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.steps-Tuple{Vector{T} where T<:Real}","page":"Home","title":"CamiXon.steps","text":"steps(x)\n\nHeatmap range transformation for steplength specification vector x\n\nExamples:\n\nx = [4,2,6]\nsteps(x)\n [0, 4, 6, 12]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.stepcenters-Tuple{Vector{T} where T<:Real}","page":"Home","title":"CamiXon.stepcenters","text":"stepcenters(x)\n\nStepcenter positions for steplength specification vector x\n\nExamples:\n\nx = [4,2,6]\nstepcenters(x)\n [2.0, 5.0, 9.0]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.stepedges-Tuple{Vector{T} where T<:Real}","page":"Home","title":"CamiXon.stepedges","text":"stepedges(x)\n\nStepedges for steplength specification vector x\n\nExamples:\n\nx = [4,2,6]\nstepedges(x)\n [0, 4, 6, 12]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.edges","page":"Home","title":"CamiXon.edges","text":"edges(px [, Δx[, x0]])\n\nHeatmap range transformation from pixel coordinates to physical coordinates, with pixelsize Δx and offset x0, both in physical units.\n\nExamples:\n\npx = 1:5\nΔx = 2.5\nx0 = 2.5\nedges(px)\n [0.5, 1.5, 2.5, 3.5, 4.5]\n\nedges(px, Δx)\n [1.25, 3.75, 6.25, 8.75, 11.25]\n\nedges(px, Δx, x0)\n [-1.25, 1.25, 3.75, 6.25, 8.75]\n\n\n\n\n\n","category":"function"},{"location":"#Search-algorithms","page":"Home","title":"Search algorithms","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"find_all(A::Union{String,AbstractArray{T,1}}, a::T...; count=false)  where T\nfind_first(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T\nfind_last(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T","category":"page"},{"location":"#CamiXon.find_all-Union{Tuple{T}, Tuple{Union{AbstractVector{T}, String}, Vararg{T, N} where N}} where T","page":"Home","title":"CamiXon.find_all","text":"find_all(A [,a...]; count=false)\n\nA: string/array of elements of the same type\n\ndefault   : Array containing the index (indices) of selected elements of A (default: all elements)\n\ncount=true: The number of indices found for selected elements of A (default: all elements)\n\nExamples:\n\nA = [:📑,:📌,:📢,:📌,:📞]\nB = [1,2,3,2,5]\nstr = \"aβcβd\";\nfind_all(A) == find_all(B) == find_all(str)\ntrue\n\nfind_all(A,:📌)\n1-element Array{Array{Int64,1},1}:\n [2, 4]\n\nfind_all(str)\n4-element Array{Array{Int64,1},1}:\n [1]\n [2, 4]\n [3]\n [5]\n\nfind_all(A; count=true)\n4-element Array{Int64,1}:\n 1\n 2\n 1\n 1\n\nstr = \"📑📌📢📌📞\"\nfind_all(str,'📌')\n1-element Array{Array{Int64,1},1}:\n [2, 4]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.find_first-Union{Tuple{T}, Tuple{Union{AbstractVector{T}, String}, Vararg{T, N} where N}} where T","page":"Home","title":"CamiXon.find_first","text":"find_first(A [,a...]; dict=false)\n\nThe first index of selected Array element\n\nA: string/array of elements of the same type\n\ndefault  : Array containing the first index (indices) of selected elements of A (default: all elements)\n\ndict=true: Dict for the first index (indices) of selected elements of A (default: all elements)\n\nExamples:\n\nA = [:📑,:📌,:📢,:📌,:📞]\nB = [1,2,3,2,5]\nstr = \"aβcβd\";\n\nfind_first(A) == find_first(B) == find_first(str)\ntrue\n\nfind_first(A,:📌)\n1-element Array{Array{Int64,1},1}:\n 2\n\nfind_last(A,:📌; dict=true)\n1-element Array{Pair{Symbol,Int64},1}:\n :📌 => 2\n\nfind_last(A; dict=true)\n4-element Array{Pair{Symbol,Int64},1}:\n :📑 => 1\n :📌 => 2\n :📢 => 3\n :📞 => 5\n\nfind_first(str)\n4-element Array{Int64,1}:\n 1\n 2\n 3\n 5\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.find_last-Union{Tuple{T}, Tuple{Union{AbstractVector{T}, String}, Vararg{T, N} where N}} where T","page":"Home","title":"CamiXon.find_last","text":"find_last(A [,a...]; dict=false)\n\nThe last index of selected Array element\n\nA: string/array of elements of the same type\n\ndefault  : Array containing the lasst index (indices) of selected elements of A (default: all elements)\n\ndict=true: Dict for the lasst index (indices) of selected elements of A (default: all elements)\n\nExamples:\n\nA = [:📑,:📌,:📢,:📌,:📞]\nB = [1,2,3,2,5]\nstr = \"aβcβd\";\nfind_last(A) == find_first(B) == find_first(str)\ntrue\n\nfind_last(A,:📌)\n1-element Array{Array{Int64,1},1}:\n 4\n\nfind_last(A,:📌; dict=true)\n1-element Array{Pair{Symbol,Int64},1}:\n :📌 => 4\n\nfind_last(A; dict=true)\n4-element Array{Pair{Symbol,Int64},1}:\n :📑 => 1\n :📌 => 4\n :📢 => 3\n :📞 => 5\n\nfind_last(str)\n4-element Array{Int64,1}:\n 1\n 4\n 3\n 5\n\n\n\n\n\n","category":"method"},{"location":"#Math","page":"Home","title":"Math","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"canonical_partitions(n::Int, m=0; header=true, reverse=true)\ninteger_partitions(n::Int, m=0; transpose=false, count=false)\nlog10_characteristic_power(x)\nlog10_mantissa(x)\npolynom_deriv_coeffs(c,deriv=0)\npolynom(c::Vector{T}, x::T) where T<:Real","category":"page"},{"location":"#CamiXon.canonical_partitions","page":"Home","title":"CamiXon.canonical_partitions","text":"canonical_partitions(n; header=false, reverse=true)\n\nThe canonical partition in integers of the integer n\n\nheader=true : unit patition included in output\n\nExamples:\n\ncanonical_partitions(6; header=true, reverse=false)\n6-element Array{Array{Int64,1},1}:\n [6]\n [5, 1]\n [4, 2]\n [3, 3]\n [2, 2, 2]\n [1, 1, 1, 1, 1, 1]\n\ncanonical_partitions(6; header=true)\n6-element Array{Array{Int64,1},1}:\n [1, 1, 1, 1, 1, 1]\n [2, 2, 2]\n [3, 3]\n [4, 2]\n [5, 1]\n [6]\n\ncanonical_partitions(6)\n5-element Array{Array{Int64,1},1}:\n [1, 1, 1, 1, 1, 1]\n [2, 2, 2]\n [3, 3]\n [4, 2]\n [5, 1]\n\n\n\n\n\n","category":"function"},{"location":"#CamiXon.integer_partitions","page":"Home","title":"CamiXon.integer_partitions","text":"integer_partitions(n [,m]; transpose=false, count=false)\n\ndefault              : The integer partitions of n\n\ncount=true           : The number of integer partitions of n\n\ntranspose=false/true : for m>0 restricted to partitions with maximum part/length m\n\ndefinitions:\n\nThe integer partition of the positive integer n is a nonincreasing sequence of positive integers p1, p2,... pk whose sum is n.\n\nThe elements of the sequence are called the parts of the partition.\n\nExamples:\n\ninteger_partitions(7)\n15-element Array{Array{Int64,1},1}:\n [1, 1, 1, 1, 1, 1, 1]\n [2, 2, 2, 1]\n [3, 3, 1]\n [4, 3]\n [5, 2]\n [6, 1]\n [7]\n [2, 2, 1, 1, 1]\n [3, 2, 2]\n [4, 2, 1]\n [5, 1, 1]\n [2, 1, 1, 1, 1, 1]\n [3, 2, 1, 1]\n [4, 1, 1, 1]\n [3, 1, 1, 1, 1]\n\ninteger_partitions(7; count=true)\n15\n\ninteger_partitions(7,4; count=true)\n3\n\ninteger_partitions(7,4)\n3-element Array{Array{Int64,1},1}:\n [4, 3]\n [4, 2, 1]\n [4, 1, 1, 1]\n\ninteger_partitions(7,4; transpose=true)\n3-element Array{Array{Int64,1},1}:\n [2, 2, 2, 1]\n [3, 2, 1, 1]\n [4, 1, 1, 1]\n\n\n\n\n\n","category":"function"},{"location":"#CamiXon.log10_characteristic_power-Tuple{Any}","page":"Home","title":"CamiXon.log10_characteristic_power","text":"log10_characteristic_power(x)\n\ncharacteristic power-of-10 of the number x\n\nExamples:\n\nlog10_characteristic_power.([3,30,300])\n3-element Vector{Int64}:\n 0\n 1\n 2\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.log10_mantissa-Tuple{Any}","page":"Home","title":"CamiXon.log10_mantissa","text":"log10_mantissa(x)\n\nlog10 mantissa of the number x\n\nExamples:\n\nlog10_mantissa.([3,30,300])\n3-element Vector{Float64}:\n 0.47712125471966244\n 0.4771212547196624\n 0.4771212547196626\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.polynom_deriv_coeffs","page":"Home","title":"CamiXon.polynom_deriv_coeffs","text":"polynom_deriv_coeffs(c[,deriv=0])\n\nCoefficients for the derivatives of the polynomial of degree d = length(c)-1 defined by the elements of the Array c[1:d+1]: p(cx) = c1 + c2 x +  + cd+1 xᵈ\n\nExamples:\n\nd = 5\nc = [1.0 for i=1:d+1]\npolynom_deriv_coeffs(c)                # default is direct copy of `c`\n [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]\npolynom_deriv_coeffs(c,1)              # coefficients of 1st derivative of `polynom(c,x)`\n [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]\npolynom_deriv_coeffs(c,2)              # coefficients of 2nd derivative of `polynom(c,x)`\n [-0.0, 0.0, 2.0, 6.0, 12.0, 20.0]\npolynom_deriv_coeffs(c,5)              # coefficients of 5th derivative of `polynom(c,x)`\n [0.0, -0.0, 0.0, -0.0, 0.0, 120.0]\npolynom_deriv_coeffs(c,6)              # coefficients of 6th derivative of `polynom(c,x)`\n [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]\n\n\n\n\n\n","category":"function"},{"location":"#CamiXon.polynom-Union{Tuple{T}, Tuple{Vector{T}, T}} where T<:Real","page":"Home","title":"CamiXon.polynom","text":"polynom(c,x)\n\nPolynomial of degree d = length(c)-1 defined by the elements of array c[1:d+1]:\n\np(cx) = c1 + c2 x +  + cd+1 xᵈ\n\nExamples:\n\nd = 5\nc = [1.0 for i=1:d+1]\nc = polynom_deriv_coeffs(c)          # default is simple\nprintln(c)\n [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]\nf(x) = polynom(c,x)\nprintln([f(1.0),f(2.0)])             # values of polynomial for x = 1.0 and x = 2.0\n [6.0, 63.0]\n\n\n\n\n\n","category":"method"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
