function test_fits_create()
    
    strExample = "minimal.fits"
    f = fits_create(strExample, protect=false)
    
    test1 = f[1].header.keys[1]  == "SIMPLE"
    test2 = f[1].dataobject.data == Any[]
    test3 = get(Dict(f[1].header.dict),"SIMPLE",0)
    test4 = get(Dict(f[1].header.dict),"NAXIS",0) == 0
    
    rm(strExample); f = nothing
    
    test = .![test1, test2, test3, test4]

    return !convert(Bool,sum(test))
    
end

function test_fits_extend()
    
    strExample = "test_example.fits"
    data = [0x0000043e, 0x0000040c, 0x0000041f]
    f = fits_create(strExample, data, protect=false)
    
    table1 = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]
    table2 = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]
    table3 = [1.23,2.12,3.,4.,5.]
    table4 = ['a','b','c','d','e']
    table5 = ["a","bb","ccc","dddd","ABCeeaeeEEEEEEEEEEEE"]
    data = [table1,table2,table3,table4,table5]
    f = fits_extend(strExample, data, "TABLE")

    test1 = f[1].header.keys[1]  == "SIMPLE"    
    test2 = f[1].dataobject.data[1] == 0x0000043e   
    test3 = f[2].header.keys[1]  == "XTENSION"  
    test4 = f[2].dataobject.data[1] == "1.0e-6 1086 1.23 a a                    "   
    test5 = get(Dict(f[2].header.dict),"NAXIS",0) == 2

    rm(strExample); f = nothing; data = nothing
    
    test = .![test1, test2, test3, test4, test5]

    return !convert(Bool,sum(test))
    
end

function test_fits_read()
    
    strExample = "minimal.fits"
    f = fits_create(strExample, protect=false)    
    f = nothing 
    f = fits_read(strExample)
    
    test1 = f[1].header.keys[1]  == "SIMPLE"
    test2 = f[1].dataobject.data == Any[]
    test3 = get(Dict(f[1].header.dict),"SIMPLE",0)
    test4 = get(Dict(f[1].header.dict),"NAXIS",0) == 0
    
    rm(strExample); f = nothing
    
    test = .![test1, test2, test3, test4]

    return !convert(Bool,sum(test))
    
end

function test_fits_info()
    
    strExample = "minimal.fits"
    f = fits_create(strExample, protect=false) 
    
    info = [
            "\r\nFile: minimal.fits",
            "hdu: 1",
            "hdutype: PRIMARY",
            "DataType: Any",
            "Datasize: (0,)",
            "\r\nMetainformation:",
            "SIMPLE  =                    T / file does conform to FITS standard             ",
            "NAXIS   =                    0 / number of data axes                            ",
            "EXTEND  =                    T / FITS dataset may contain extensions            ",
            "COMMENT    Primary FITS HDU    / http://fits.gsfc.nasa.gov/iaufwg               ",
            "END                                                                             "
            ]
       
    test = fits_info(f[1]; printformat=false) == Base.join(info .* "\r\n")
    
    return test
    
end
