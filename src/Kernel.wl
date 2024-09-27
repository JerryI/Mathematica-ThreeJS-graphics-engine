BeginPackage["JerryI`Notebook`Graphics3DUtils`", {
  "JerryI`Misc`Events`", "Notebook`Editor`Kernel`FrontSubmitService`",
  "Notebook`Editor`FrontendObject`", 
  "Notebook`Editor`Boxes`"
}]

Metalness::usage = "Specify metallness of the surface Metalness[1] used in Graphics3D"
Emissive::usage = "Makes a surface emitt light Emissive[RGBColor[...], intensity_:1] used in Graphics3D"
Roughness::usage = "Specify the roughness of the surface Roughness[1] used in Graphics3D"
Shadows::usage = "used in Graphics3D. Decide if you need to cast shadows from objects. Shadows[True]"

HemisphereLight::usage = "HemisphereLight[skyColor_RGBColor, groundColor_RGBColor, intensity_] used in Graphics3D"

MeshMaterial::usage = "specifies the material for 3D primitives. MeshMaterial[MeshPhysicalMaterial[]], MeshMaterial[MeshToonMaterial[]]"

MeshPhysicalMaterial::usage = ""
MeshToonMaterial::usage = ""

EventListener::usage = "Internal wrapper for Graphics object to catch events"

Graphics3D`Materials;
Graphics3D`Serialize;

Begin["`Private`"]

Graphics3D`Materials["Glass"] = Directive[White, "EmissiveIntensity"->0, "Ior"->1.51, "Transmission"->1.0, "Roughness"->0.13]
Graphics3D`Materials["Iridescent"] = Directive[RGBColor["#474747"], "Roughness"->0.25, "Metalness"->1.0, "Iridescence"->1.0, "IridescenceIOR"->2.2]
Graphics3D`Materials["Acrylic"] = Directive[White, "Roughness"->0, "Metalness"->0, "Transmission"->1.0, , "AttenuationDistance"->0.75, "AttenuationColor"->RGBColor["#2a6dc6"] ]

(*
    list of all properties supported in Directives

    "Color",
    "Emissive",
    "Emissiveintensity",
    "Roughness",
    "Metalness",
    "Ior",
    "Transmission",
    "Thinfilm",
    "MaterialThickness",
    "Attenuationcolor",
    "Attenuationdistance",
    "Opacity",
    "Clearcoat",
    "Clearcoatroughness",
    "Sheencolor",
    "Sheenroughness",
    "Iridescence",
    "Iridescenceior",
    "Iridescencethickness",
    "Specularcolor",
    "Specularintensity",
    "Matte",
    "Flatshading",
    "Castshadow"

*)

listener[p_, list_] := With[{uid = CreateUUID[]}, With[{
    rules = Map[Function[rule, rule[[1]] -> uid ], list]
},
    EventHandler[uid, list];
    EventListener[p, rules]
] ]

Unprotect[Sphere];

Sphere      /: EventHandler[p_Sphere, list_List] := listener[p, list]

Protect[Sphere];


(* CALL A MODAL HERE 
Unprotect[Rasterize]
Rasterize[g_Graphics3D, any___] := With[{base = FrontFetch[Graphics3D`Serialize[Plot3D[Sin[x + y^2], {x, -3, 3}, {y, -2, 2}], "TemporalDOM"->True] ]},
  ImportString[StringDrop[base, StringLength["data:image/png;base64,"] ], "Base64"]
]
*)

Unprotect[Graphics3D]

System`WLXForm;

Graphics3D /: MakeBoxes[g_Graphics3D, WLXForm] := With[{fe = CreateFrontEndObject[g]},
    MakeBoxes[fe, WLXForm]
]

Graphics3D /: MakeBoxes[System`Dump`g_Graphics3D,System`Dump`fmt:StandardForm|TraditionalForm] := If[ByteCount[System`Dump`g] < 2 1024,
    ViewBox[System`Dump`g, System`Dump`g]
,
    With[{fe = CreateFrontEndObject[System`Dump`g]},
        MakeBoxes[fe, System`Dump`fmt]
    ]
]

End[]
EndPackage[]