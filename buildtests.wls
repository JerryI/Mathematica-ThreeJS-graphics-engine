Import["https://raw.githubusercontent.com/JerryI/tinyweb-mathematica/master/WSP/WSP.wl"];

JerryI`WolframJSFrontend`root = FileNameJoin[{"tests"}];
Get[FileNameJoin[{"tests", "modules", "Packages.wl"}]];

LoadPluginsConfiguration;

ExportFunction[func_] := Module[{name = FileBaseName[filename], json = ExportString[func, "ExpressionJSON"]},
    Export["tests/public/js/"<>name<>".js", StringTemplate["var `` = ``;"][name, json], "String"];
    Export["tests/public/imgs/"<>name<>".png", Show[func, ImageSize->400], "PNG"];
    json
];

$publicpath = "tests";

With[{name = FileBaseName[#]},

    Export["tests/public/"<>name<>".html", LoadPage["src/"<>FileNameTake[#], {filename = #, function = Null, json = ""}], "String"];
    StringTemplate["exported ``.html"][name]//Print;

] &/@ SortBy[FileNames["*.wsp", "tests/src"], StringTake[FileBaseName[#], 1] == "_" &];