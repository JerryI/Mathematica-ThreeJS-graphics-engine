Import["https://raw.githubusercontent.com/JerryI/tinyweb-mathematica/master/WSP/WSP.wl"];

ExportFunction[func_] := Module[{name = FileBaseName[filename], json = ExportString[func, "ExpressionJSON"]},
    Export["tests/js/"<>name<>".js", StringTemplate["var `` = ``;"][name, json], "String"];
    Export["tests/imgs/"<>name<>".png", Show[func, ImageSize->400], "PNG"];
    json
];

$publicpath = "tests";

With[{name = FileBaseName[#]},

    Export["tests/"<>name<>".html", LoadPage["src/"<>FileNameTake[#], {filename = #, function = Null, json = ""}], "String"];
    StringTemplate["exported ``.html"][name]//Print;

] &/@ SortBy[FileNames["*.wsp", "tests/src"], StringTake[FileBaseName[#], 1] == "_" &];