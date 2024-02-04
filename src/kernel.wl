BeginPackage["JerryI`Notebook`Graphics3DUtils`", {"JerryI`Misc`Events`"}]

Metalness::usage = "Specify metallness of the surface Metalness[1] used in Graphics3D"
Emissive::usage = "Makes a surface emitt light Emissive[RGBColor[...], intensity_:1] used in Graphics3D"
Roughness::usage = "Specify the roughness of the surface Roughness[1] used in Graphics3D"
Shadows::usage = "used in Graphics3D. Decide if you need to cast shadows from objects. Shadows[True]"

EventListener::usage = "Internal wrapper for Graphics object to catch events"

Begin["`Private`"]

listener[p_, list_] := With[{uid = CreateUUID[]}, With[{
    rules = Map[Function[rule, rule[[1]] -> uid ], list]
},
    EventHandler[uid, list];
    EventListener[p, rules]
] ]

Unprotect[Sphere];

Sphere      /: EventHandler[p_Point, list_List] := listener[p, list]

Protect[Sphere];

End[]
EndPackage[]