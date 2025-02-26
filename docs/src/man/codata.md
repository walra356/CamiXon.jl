## Codata

```@docs
Codata
```
#### castCodata

```@docs
castCodata(year::Int)
```

#### listCodata

```@docs
listCodata(codata::Codata; msg=true)
```

## Value

```@docs
Value
```

#### strValue

```@docs
strValue(f::Value)
```

## NamedValue

```@docs
NamedValue
```

#### castNamedValue

```@docs
castNamedValue(val::Value; name=" ", comment=" ")
```

## Unit conversion

```@docs
convertUnit(val, codata; unitIn="Hartree", unitOut="xHz")
```

## Calibration report

```@docs
calibrationReport(E, Ecal, codata::Codata; unitIn="Hartree", msg=true)
```
