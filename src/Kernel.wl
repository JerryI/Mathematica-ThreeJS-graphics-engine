BeginPackage["JerryI`Notebook`Graphics3DUtils`", {"JerryI`Misc`Events`"}]

Metalness::usage = "Specify metallness of the surface Metalness[1] used in Graphics3D"
Emissive::usage = "Makes a surface emitt light Emissive[RGBColor[...], intensity_:1] used in Graphics3D"
Roughness::usage = "Specify the roughness of the surface Roughness[1] used in Graphics3D"
Shadows::usage = "used in Graphics3D. Decide if you need to cast shadows from objects. Shadows[True]"

HemisphereLight::usage = "HemisphereLight[skyColor_RGBColor, groundColor_RGBColor, intensity_] used in Graphics3D"

MeshMaterial::usage = "specifies the material for 3D primitives. MeshMaterial[MeshPhysicalMaterial[]], MeshMaterial[MeshToonMaterial[]]"

MeshPhysicalMaterial::usage = ""
MeshToonMaterial::usage = ""

EventListener::usage = "Internal wrapper for Graphics object to catch events"

Graphics3D`Serialize;

Begin["`Private`"]

listener[p_, list_] := With[{uid = CreateUUID[]}, With[{
    rules = Map[Function[rule, rule[[1]] -> uid ], list]
},
    EventHandler[uid, list];
    EventListener[p, rules]
] ]

Unprotect[Sphere];

Sphere      /: EventHandler[p_Sphere, list_List] := listener[p, list]

Protect[Sphere];


(*

Then[FrontFetchAsync[Graphics3D`Serialize[Plot3D[Sin[x + y^2], {x, -3, 3}, {y, -2, 2}], "TemporalDOM"->True]], Function[val,
  test = val;
]]

ImportString[
 StringDrop[test, StringLength["data:image/png;base64,"]], 
 "Base64"]

*)

End[]
EndPackage[]