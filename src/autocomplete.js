window.EditorAutocomplete.extend([  
    {
        "label": "PointerLockControls",
        "type": "keyword",
        "info": '"Controls"->"PointerLockControls" is an option for Graphics3D to use FPS-like controls'  
    },  
    {
        "label": "Controls",
        "type": "keyword",
        "info": "Options for controls of camera in Graphics3D"  
    },
    {
        "label": "RTX",
        "type": "keyword",
        "info": '"RTX"->True is an option for Graphics3D to switch to Path-tracing engine'  
    },
    {
        "label": 'Lightmap',
        "type": 'keyword',
        'info': 'An option for Graphics3D to use HDRI light maps. Specify the url Lightmap->"https://..."'
    },
    {
        "label": "PointLight",
        "type": "keyword",
        "info": 'A point light object used in Graphics3D. Default parameters are PointLight[color=White, intensity=1, distance=0, decay=2, "Position"->{0,0,10}]'  
    },
    {
        "label": "HemisphereLight",
        "type": "keyword",
        "info": 'A hemisphere light object used in Graphics3D. Default parameters are HemisphereLight[skyColor=0xffffbb, groundColor=0x080820, intensity=1]'  
    },
    {
        "label": "MeshMaterial",
        "type": "keyword",
        "info": 'Specify mesh material used in Graphics3D for lighting. Typical parameters: MeshPhysicalMaterial, MeshLambertMaterial, MeshPhongMaterial, MeshToonMaterial'  
    },
    {
        "label": "MeshPhysicalMaterial",
        "type": "keyword",
        "info": 'A mesh material used in MeshMaterial[]'  
    },
    {
        "label": "MeshLambertMaterial",
        "type": "keyword",
        "info": 'A mesh material used in MeshMaterial[]'  
    },
    {
        "label": "MeshPhongMaterial",
        "type": "keyword",
        "info": 'A mesh material used in MeshMaterial[]'  
    },
    {
        "label": "MeshToonMaterial",
        "type": "keyword",
        "info": 'A mesh material used in MeshMaterial[]'  
    },
    {
        "label": "Metalness",
        "type": "keyword",
        "info": 'Specify metallness of the surface Metalness[1] used in Graphics3D'  
    },
    {
        "label": "Emissive",
        "type": "keyword",
        "info": 'Makes a surface emitt light Emissive[RGBColor[...]] used in Graphics3D'  
    },
    {
        "label": "Roughness",
        "type": "keyword",
        "info": 'Specify the roughness of the surface Roughness[1] used in Graphics3D'  
    },
    {
        "label": "Reflectivity",
        "type": "keyword",
        "info": 'Specify the reflectivity of the surface Reflectivity[1] used in Graphics3D'  
    },
    {
        "label": "Clearcoat",
        "type": "keyword",
        "info": 'Subsurface scattering Clearcoat[1] used in Graphics3D'  
    }                              
])

console.log('loaded!');