	const flatten = (ary) => ary.reduce((a, b) => a.concat(Array.isArray(b) ? flatten(b) : b), [])
	
	function computeGroupCenter(group) {
		var center = new THREE.Vector3();
		var children = group.children;
		var count = children.length;
		for (var i = 0; i < count; i++) {
			center.add(children[i].position);
		}
		center.divideScalar(count);
		return center;
	}
	
	function RGBtoColor(i,k,j) {
		var r = Math.round(255*i);
		var g = Math.round(255*k);
		var b = Math.round(255*j);
		
		return(new THREE.Color("rgb("+r+","+g+","+b+")"));
		
	}

	var interpretate = function (d, parent=undefined, params=undefined, mesh=undefined) {
		
		if (typeof d !== 'object') return d;
		
		var func = {
			name: d[0],
			args: d.slice(1,d.length),
			parent: parent
		};
	
		
		switch(func.name) {
			case 'Graphics3D':
				console.log("GRAPHICS3D");
				
				var data = {"axes": {}, 
							"extent": {"zmax": 1.0, "ymax": 1.0, "zmin": -1.0, "xmax": 1.0, "xmin": -1.0, "ymin": -1.0}, 

							"lighting": [ {"type": "Ambient", "color": [0.3, 0.2, 0.4]},
											{"type": "Directional", "color": [0.8, 0., 0.],
											"position": [2, 0, 2]},
											{"type": "Directional", "color": [0., 0.8, 0.],
											"position": [2, 2, 2]},
											{"type": "Directional", "color": [0., 0., 0.8],
											"position": [0, 2, 2]}],
				 
							"viewpoint":[1.3, -2.4, 2]};
							
				var container = document.getElementById("graphics3d");
				
                var camera, scene, renderer, boundbox, hasaxes, viewpoint,
                  isMouseDown = false, onMouseDownPosition,
                  tmpx, tmpy, tmpz, 
                  theta, onMouseDownTheta, phi, onMouseDownPhi;

                // Scene
                scene = new THREE.Scene();
				
				var group = new THREE.Group();
				
				params = {
	
					matrix: new THREE.Matrix4(),
					color: RGBtoColor(1,1,1),
					opacity: 1,
					thickness: 1,
					edgecolor: RGBtoColor(0,0,0)
				};
				
				params.matrix.set( 1,0,0,0,
								0,1,0,0,
								0,0,1,0,
								0,0,0,1 );
				
				interpretate(func.args[0], func, params, group);
				
				var bbox = new THREE.Box3().setFromObject(group);
					
				var center = new THREE.Vector3().addVectors(bbox.max,bbox.min).divideScalar(2);
				console.log("BBOX CENTER");
				console.log(center);
				console.log(bbox);
				//var	translate = new THREE.Matrix4().makeTranslation(-center.x,-center.y,-center.z,1);
				//group.applyMatrix(translate);	
				scene.position = center;
				
								
				  
                // Center of the scene
                //var center = new THREE.Vector3(
                //  0.5 * (data.extent.xmin + data.extent.xmax),
                //  0.5 * (data.extent.ymin + data.extent.ymax), 
                //  0.5 * (data.extent.zmin + data.extent.zmax));
                
                // Where the camera is looking
                var focus = new THREE.Vector3(center.x, center.y, center.z);
                
                // Viewpoint
                viewpoint = new THREE.Vector3(data.viewpoint[0], data.viewpoint[1], data.viewpoint[2]).sub(focus);
                
				var ln = new THREE.Vector3().addVectors(bbox.max,bbox.min.clone().negate()).length();
				
				console.log("Radius is "+ln);
				
				viewpoint.x *= ln; 
				viewpoint.y *= ln; 
				viewpoint.z *= ln;
				
				var radius = viewpoint.length()
                
                onMouseDownTheta = theta = Math.acos(viewpoint.z / radius);
                onMouseDownPhi = phi = (Math.atan2(viewpoint.y, viewpoint.x) + 2*Math.PI) % (2 * Math.PI);
           
                
                
                
                camera = new THREE.PerspectiveCamera(
                  35,             // Field of view
                  1,            // Aspect ratio
                  0.1*radius,     // Near plane
                  1000*radius     // Far plane
                );
                
                function update_camera_position() {
                  camera.position.x = radius * Math.sin(theta) * Math.cos(phi);
                  camera.position.y = radius * Math.sin(theta) * Math.sin(phi);
                  camera.position.z = radius * Math.cos(theta);
                  camera.position.add(focus);
                  camera.lookAt(focus);
                }
                
                update_camera_position();
                camera.up = new THREE.Vector3(0,0,1);
                
                scene.add(camera);
                
                // Lighting
                function addLight(l) {
                  var color = new THREE.Color().setRGB(l.color[0], l.color[1], l.color[2]);
                  var light;
                
                  if (l.type === "Ambient") {
                    light = new THREE.AmbientLight(color.getHex());
                  } else if (l.type === "Directional") {
                    light = new THREE.DirectionalLight(color.getHex(), 1);
                  } else if (l.type === "Spot") {
                    light = new THREE.SpotLight(color.getHex());
                    light.position.set(l.position[0], l.position[1], l.position[2]);
                    light.target.position.set(l.target[0], l.target[1], l.target[2]);
                    light.target.updateMatrixWorld(); // This fixes bug in THREE.js
                    light.angle = l.angle;
                  } else if (l.type === "Point") {
                    light = new THREE.PointLight(color.getHex());
                    light.position.set(l.position[0], l.position[1], l.position[2]);
                
                    // Add visible light sphere
                    var lightsphere = new THREE.Mesh(
                      new THREE.SphereGeometry(0.007*radius, 16, 8),
                      new THREE.MeshBasicMaterial({color: color.getHex()})
                    );
                    lightsphere.position = light.position;
                    scene.add(lightsphere);
                  } else {
                    alert("Error: Internal Light Error", l.type);
                    return;
                  }
                  return light;
                }
                
                function getInitLightPos(l) {
                  // Initial Light position in spherical polar coordinates
                  if (l.position instanceof Array) {
                    var tmppos = new THREE.Vector3(l.position[0], l.position[1], l.position[2]);
                    var result = {"radius": radius * tmppos.length()};
                
                    if (tmppos.length() <= 0.0001) {
                      result.theta = 0;
                      result.phi = 0;
                    } else {
                      result.phi = (Math.atan2(tmppos.y, tmppos.x) + 2 * Math.PI) % (2 * Math.PI);
                      result.theta = Math.asin(tmppos.z / result.radius);
                    }
                    return result;
                  }
                  return;
                }
                
                function positionLights() {
                  for (var i = 0; i < lights.length; i++) {
                    if (lights[i] instanceof THREE.DirectionalLight) {
                      lights[i].position.x = initLightPos[i].radius * Math.sin(theta + initLightPos[i].theta) * Math.cos(phi + initLightPos[i].phi);
                      lights[i].position.y = initLightPos[i].radius * Math.sin(theta + initLightPos[i].theta) * Math.sin(phi + initLightPos[i].phi);
                      lights[i].position.z = initLightPos[i].radius * Math.cos(theta + initLightPos[i].theta);
                      lights[i].position.add(focus);
                    }
                  }
                }
                
                var lights = new Array(data.lighting.length);
                var initLightPos = new Array(data.lighting.length);
                
                for (var i = 0; i < data.lighting.length; i++) {
                  initLightPos[i] = getInitLightPos(data.lighting[i]);
                  
                  lights[i] = addLight(data.lighting[i]);
                  scene.add(lights[i]);
                }
                
                // BoundingBox
                boundbox = new THREE.Mesh(
                  new THREE.CubeGeometry(
					bbox.max.x - bbox.min.x,
					bbox.max.y - bbox.min.y,
					bbox.max.z - bbox.min.z),
                  new THREE.MeshBasicMaterial({color: 0x666666, wireframe: true})
                );
                boundbox.position = center;

				var geo = new THREE.EdgesGeometry( new THREE.CubeGeometry(
					bbox.max.x - bbox.min.x,
					bbox.max.y - bbox.min.y,
					bbox.max.z - bbox.min.z) ); // or WireframeGeometry( geometry )
				
				var mat = new THREE.LineBasicMaterial( { color: 0x666666, linewidth: 2 } );
				
				var wireframe = new THREE.LineSegments( geo.translate(center.x,center.y,center.z), mat );				
				
				
                //scene.add(wireframe);  
                
                // Draw the Axes
                if (data.axes.hasaxes instanceof Array) {
                  hasaxes = new Array(data.axes.hasaxes[0], data.axes.hasaxes[1], data.axes.hasaxes[2]);
                } else if (data.axes.hasaxes instanceof Boolean) {
                  if (data.axes) {
                    hasaxes = new Array(true, true, true);
                  } else {
                    hasaxes = new Array(false, false, false);
                  }
                } else {
                  hasaxes = new Array(false, false, false);
                }
                var axesmat = new THREE.LineBasicMaterial({ color: 0x000000, linewidth : 1.5 });
                var axesgeom = [];
                var axesindicies = [
                  [[0,5], [1,4], [2,7], [3,6]],
                  [[0,2], [1,3], [4,6], [5,7]],
                  [[0,1], [2,3], [4,5], [6,7]]
                ];
                
                var axesmesh = new Array(3);
                for (var i=0; i<3; i++) {
                  if (hasaxes[i]) {
                    axesgeom[i] = new THREE.Geometry();
                    axesgeom[i].vertices.push(new THREE.Vector3().addVectors(
                      boundbox.geometry.vertices[axesindicies[i][0][0]], boundbox.position)
                    );
                    axesgeom[i].vertices.push(new THREE.Vector3().addVectors(
                      boundbox.geometry.vertices[axesindicies[i][0][1]], boundbox.position)
                    );
                    axesmesh[i] = new THREE.Line(axesgeom[i], axesmat);
                    axesmesh[i].geometry.dynamic = true;
                    scene.add(axesmesh[i]);
                  }
                }
                
                function boxEdgeLength(i, j) {
                  edge = new THREE.Vector3().sub(
                    toCanvasCoords(boundbox.geometry.vertices[axesindicies[i][j][0]]),
                    toCanvasCoords(boundbox.geometry.vertices[axesindicies[i][j][1]])
                  );
                  edge.z = 0;
                  return edge.length();
                }
                
                function positionAxes() {
                  // Automatic axes placement
                  nearj = null;
                  nearl = 10*radius;
                  farj = null;
                  farl = 0.0;
                  
                  tmpv = new THREE.Vector3();
                  for (var j = 0; j < 8; j++) {
                    tmpv.addVectors(boundbox.geometry.vertices[j], boundbox.position);
                    tmpv.sub(camera.position);
                    tmpl = tmpv.length();
                    if (tmpl < nearl) {
                      nearl = tmpl;
                      nearj = j;
                    } else if (tmpl > farl) {
                      farl = tmpl;
                      farj = j;
                    }
                  }
                  for (var i = 0; i < 3; i++) {
                    if (hasaxes[i]) {
                      maxj = null;
                      maxl = 0.0;
                      for (var j = 0; j < 4; j++) {
                        if (axesindicies[i][j][0] !== nearj && axesindicies[i][j][1] !== nearj && axesindicies[i][j][0] !== farj && axesindicies[i][j][1] !== farj) {
                          tmpl = boxEdgeLength(i, j);
                          if (tmpl > maxl) {
                            maxl = tmpl;
                            maxj = j;
                          }
                        }
                      }
                      axesmesh[i].geometry.vertices[0].addVectors(boundbox.geometry.vertices[axesindicies[i][maxj][0]], boundbox.position);
                      axesmesh[i].geometry.vertices[1].addVectors(boundbox.geometry.vertices[axesindicies[i][maxj][1]], boundbox.position);
                      axesmesh[i].geometry.verticesNeedUpdate = true;
                    }
                  }
                  update_axes();
                }
                
                // Axes Ticks
                var tickmat = new THREE.LineBasicMaterial({ color: 0x000000, linewidth : 1.2 });
                var ticks = new Array(3);
                var ticks_small = new Array(3);
                var ticklength = 0.005*radius;
                
                for (var i = 0; i < 3; i++) {
                  if (hasaxes[i]) {
                    ticks[i] = [];
                    for (var j = 0; j < data.axes.ticks[i][0].length; j++) {
                      tickgeom = new THREE.Geometry();
                      tickgeom.vertices.push(new THREE.Vector3());
                      tickgeom.vertices.push(new THREE.Vector3());
                      ticks[i].push(new THREE.Line(tickgeom, tickmat));
                      scene.add(ticks[i][j]);
                
                    }
                    ticks_small[i] = [];
                    for (var j = 0; j < data.axes.ticks[i][1].length; j++) {
                       tickgeom = new THREE.Geometry();
                       tickgeom.vertices.push(new THREE.Vector3());
                       tickgeom.vertices.push(new THREE.Vector3());
                       ticks_small[i].push(new THREE.Line(tickgeom, tickmat));
                       scene.add(ticks_small[i][j]);
                    }
                  }
                }
                
                function getTickDir(i) {
                  var tickdir = new THREE.Vector3();
                  if (i === 0) {
                    if (0.25*Math.PI < theta && theta < 0.75*Math.PI) {
                      if (axesgeom[0].vertices[0].z > boundbox.position.z) {
                        tickdir.set(0, 0, -ticklength);
                      } else {
                        tickdir.set(0, 0, ticklength);
                      }
                    } else {
                      if (axesgeom[0].vertices[0].y > boundbox.position.y) {
                        tickdir.set(0,-ticklength, 0);
                      } else {
                        tickdir.set(0, ticklength, 0);
                      }
                    }
                  } else if (i === 1) {
                    if (0.25*Math.PI < theta && theta < 0.75*Math.PI) {
                      if (axesgeom[1].vertices[0].z > boundbox.position.z) {
                        tickdir.set(0, 0, -ticklength);
                      } else {
                        tickdir.set(0, 0, ticklength);
                      }
                    } else {
                      if (axesgeom[1].vertices[0].x > boundbox.position.x) {
                        tickdir.set(-ticklength, 0, 0);
                      } else {
                        tickdir.set(ticklength, 0, 0);
                      }
                    }
                  } else if (i === 2) {
                    if ((0.25*Math.PI < phi && phi < 0.75*Math.PI) || (1.25*Math.PI < phi && phi < 1.75*Math.PI)) {
                      if (axesgeom[2].vertices[0].x > boundbox.position.x) {
                        tickdir.set(-ticklength, 0, 0);
                      } else {
                        tickdir.set(ticklength, 0, 0);
                      }
                    } else {
                      if (axesgeom[2].vertices[0].y > boundbox.position.y) {
                        tickdir.set(0, -ticklength, 0, 0);
                      } else {
                        tickdir.set(0, ticklength, 0, 0);
                      }
                    }
                  }
                  return tickdir;
                }
                
                function update_axes() {
                  for (var i = 0; i < 3; i++) {
                    if (hasaxes[i]) {
                      tickdir = getTickDir(i);
                      small_tickdir = tickdir.clone();
                      small_tickdir.multiplyScalar(0.5);
                      for (var j = 0; j < data.axes.ticks[i][0].length; j++) {
                        tmpval = data.axes.ticks[i][0][j];
                
                        ticks[i][j].geometry.vertices[0].copy(axesgeom[i].vertices[0]);
                        ticks[i][j].geometry.vertices[1].addVectors(axesgeom[i].vertices[0], tickdir);
                
                        if (i === 0) {
                          ticks[i][j].geometry.vertices[0].x = tmpval;
                          ticks[i][j].geometry.vertices[1].x = tmpval;
                        } else if (i === 1) {
                          ticks[i][j].geometry.vertices[0].y = tmpval;
                          ticks[i][j].geometry.vertices[1].y = tmpval;
                        } else if (i === 2) {
                          ticks[i][j].geometry.vertices[0].z = tmpval;
                          ticks[i][j].geometry.vertices[1].z = tmpval;
                        }
                
                        ticks[i][j].geometry.verticesNeedUpdate = true;
                      }
                      for (var j = 0; j < data.axes.ticks[i][1].length; j++) {
                        tmpval = data.axes.ticks[i][1][j];
                
                        ticks_small[i][j].geometry.vertices[0].copy(axesgeom[i].vertices[0]);
                        ticks_small[i][j].geometry.vertices[1].addVectors(axesgeom[i].vertices[0], small_tickdir);
                
                        if (i === 0) {
                          ticks_small[i][j].geometry.vertices[0].x = tmpval;
                          ticks_small[i][j].geometry.vertices[1].x = tmpval;
                        } else if (i === 1) {
                          ticks_small[i][j].geometry.vertices[0].y = tmpval;
                          ticks_small[i][j].geometry.vertices[1].y = tmpval;
                        } else if (i === 2) {
                          ticks_small[i][j].geometry.vertices[0].z = tmpval;
                          ticks_small[i][j].geometry.vertices[1].z = tmpval;
                        }
                
                        ticks_small[i][j].geometry.verticesNeedUpdate = true;
                      }
                    }
                  }
                }
                update_axes();
                
                // Axes numbering using divs
                var ticknums = new Array(3);
                for (var i = 0; i < 3; i++) {
                  if (hasaxes[i]) {
                    ticknums[i] = new Array(data.axes.ticks[i][0].length);
                    for (var j = 0; j < ticknums[i].length; j++) {
                      ticknums[i][j] = document.createElement('div');
                      ticknums[i][j].innerHTML = data.axes.ticks[i][2][j];
                
                      // Handle Minus signs
                      if (data.axes.ticks[i][0][j] >= 0) {
                        ticknums[i][j].style.paddingLeft = "0.5em";
                      } else {
                        ticknums[i][j].style.paddingLeft = 0;
                      }
                
                      ticknums[i][j].style.position = "absolute";
                      ticknums[i][j].style.fontSize = "0.8em";
                      container.appendChild(ticknums[i][j]);
                    }
                  }
                }
                
                function toCanvasCoords(position) {
                  var pos = position.clone();
                  var projScreenMat = new THREE.Matrix4();
                  projScreenMat.multiply(camera.projectionMatrix, camera.matrixWorldInverse);
                  //.multiplyVector3( pos );
				  pos = pos.applyMatrix4(projScreenMat);
                
                  var result = new THREE.Vector3((pos.x + 1 ) * 200, (1-pos.y) * 200, (pos.z + 1 ) * 200);
                  return result;
                }
                
                function positionticknums() {
                  for (var i = 0; i < 3; i++) {
                    if (hasaxes[i]) {
                      for (var j = 0; j < ticknums[i].length; j++) {
                        var tickpos3D = ticks[i][j].geometry.vertices[0].clone();
                        var tickDir = new THREE.Vector3().sub(ticks[i][j].geometry.vertices[0], ticks[i][j].geometry.vertices[1]);
                        //tickDir.multiplyScalar(3);
                        tickDir.setLength(3*ticklength);
                        tickDir.x *= 2.0;
                        tickDir.y *= 2.0;
                        tickpos3D.add(tickDir);
                        var tickpos = toCanvasCoords(tickpos3D);
                        tickpos.x -= 10;
                        tickpos.y += 8;
                
                        ticknums[i][j].style.left = tickpos.x.toString() + "px";
                        ticknums[i][j].style.top = tickpos.y.toString() + "px";
                        if (tickpos.x < 5 || tickpos.x > 395 || tickpos.y < 5 || tickpos.y > 395) {
                          ticknums[i][j].style.display = "none";
                        }
                        else {
                          ticknums[i][j].style.display = "";
                        }
                      }
                    }
                  }
                }
                
                
			
                scene.add(group);
				//loader = new THREE.JSONLoader();

				//loader.load( JSON.parse(str));
				//loader.onLoadComplete=function(mesh){scene.add( mesh )} 
                
                // Plot the primatives
                /*for (var indx = 0; indx < data.elements.length; indx++) {
                  var type = data.elements[indx].type;
                  switch(type) {
                    case "point":
                      scene.add(drawPoint(data.elements[indx]));
                      break;
                    case "line":
                      scene.add(drawLine(data.elements[indx]));
                      break;
                    case "polygon":
                      scene.add(drawPolygon(data.elements[indx]));
                      break;
                    case "sphere":
                      scene.add(drawSphere(data.elements[indx]));
                      break;
                    case "cube":
                      scene.add(drawCube(data.elements[indx]));
                      break;
                    default:
                      alert("Error: Unknown type passed to drawGraphics3D");
                  }
                }*/
                
                // Renderer (set preserveDrawingBuffer to deal with issue
                // of weird canvas content after switching windows)
                if (Detector.webgl) {
                  renderer = new THREE.WebGLRenderer({antialias: true, preserveDrawingBuffer: true});
                } else { 
                  renderer = new THREE.CanvasRenderer({antialias: true, preserveDrawingBuffer: true});
                
                  message = document.createElement('div');
                  message.innerHTML = "Canvas Renderer support is experimental, please enable WebGL where possible.";
                  message.style.position = "absolute";
                  message.style.fontSize = "0.8em";
                  message.style.color = "#FF6060";
                  container.appendChild(message);
                }
                
                renderer.setSize(400, 400);
                renderer.setClearColor( 0xffffff );
                container.appendChild(renderer.domElement);
                
                function render() {
                  positionLights();
                  renderer.render( scene, camera );
                }
                
                function toScreenCoords(position) {
                  return position.clone().applyMatrix3(camera.matrixWorldInverse);
                }
                
                function ScaleInView() {
                  var tmp_fov = 0.0;
                  //var proj2d = new THREE.Vector3();
                
                  /*for (var i=0; i<boundbox.geometry.vertices.length; i++) {
                    proj2d = proj2d.addVectors(boundbox.geometry.vertices[i], boundbox.position);
                    proj2d = toScreenCoords(proj2d);
                
                    angle = 114.59 * Math.max(
                       Math.abs(Math.atan(proj2d.x/proj2d.z) / camera.aspect),
                       Math.abs(Math.atan(proj2d.y/proj2d.z))
                    );
                    tmp_fov = Math.max(tmp_fov, angle);
                  }*/
				  //console.log(bbox);
				  var height = bbox.min.clone().sub(bbox.max).length();
				  var dist = center.clone().sub(camera.position).length();
				  tmp_fov = 2 * Math.atan( height / ( 2 * dist ) ) * ( 180 / Math.PI );
				  
                  camera.fov = tmp_fov + 5;
                  camera.updateProjectionMatrix();
                }
                
                // Mouse Interactions
                function onDocumentMouseDown( event ) {
                  event.preventDefault();
                
                  isMouseDown = true;
                  isShiftDown = false;
                  isCtrlDown = false;
                
                  onMouseDownTheta = theta;
                  onMouseDownPhi = phi;
                
                  onMouseDownPosition.x = event.clientX;
                  onMouseDownPosition.y = event.clientY;
                
                  onMouseDownFocus = new THREE.Vector3().copy(focus);
                }
                
                function onDocumentMouseMove(event) {
                  event.preventDefault();
                
                  if (isMouseDown) {
                    positionticknums();
                
                    if (event.shiftKey) {
                      // console.log("Pan");
                      if (! isShiftDown) {
                        isShiftDown = true;
                        onMouseDownPosition.x = event.clientX;
                        onMouseDownPosition.y = event.clientY;
                        autoRescale = false;
                        container.style.cursor = "move";
                      }
                      var camz = new THREE.Vector3().sub(focus, camera.position);
                      camz.normalize();
                
                      var camx = new THREE.Vector3(
                          - radius * Math.cos(theta) * Math.sin(phi) * (theta<0.5*Math.PI?1:-1),
                          radius * Math.cos(theta) * Math.cos(phi) * (theta<0.5*Math.PI?1:-1),
                          0
                      );
                      camx.normalize();
                
                      var camy = new THREE.Vector3();
                      camy.cross(camz, camx);
                
                      focus.x = onMouseDownFocus.x + (radius / 400)*(camx.x * (onMouseDownPosition.x - event.clientX) + camy.x * (onMouseDownPosition.y - event.clientY));
                      focus.y = onMouseDownFocus.y + (radius / 400)*(camx.y * (onMouseDownPosition.x - event.clientX) + camy.y * (onMouseDownPosition.y - event.clientY));
                      focus.z = onMouseDownFocus.z + (radius / 400)*(camx.z * (onMouseDownPosition.x - event.clientX) + camy.z * (onMouseDownPosition.y - event.clientY));
                
                      update_camera_position();
                
                    } else if (event.ctrlKey) {
                      // console.log("Zoom");
                      if (! isCtrlDown) {
                        isCtrlDown = true;
                        onCtrlDownFov = camera.fov;
                        onMouseDownPosition.x = event.clientX;
                        onMouseDownPosition.y = event.clientY;
                        autoRescale = false;
                        container.style.cursor = "crosshair";
                      }
                      camera.fov =  onCtrlDownFov + 20 * Math.atan((event.clientY - onMouseDownPosition.y)/50);
					  
                      camera.fov = Math.max(1, Math.min(camera.fov, 150));
					  //console.log("fov"+camera.fov);
                      camera.updateProjectionMatrix();
					  //console.log(JSON.stringify(camera));
                    } else {
                      // console.log("Spin");
                      if (isCtrlDown || isShiftDown) {
                        onMouseDownPosition.x = event.clientX;
                        onMouseDownPosition.y = event.clientY;
                        isShiftDown = false;
                        isCtrlDown = false;
                        container.style.cursor = "pointer";
                      }
                
                      phi = 2 * Math.PI * (onMouseDownPosition.x - event.clientX) / 400 + onMouseDownPhi;
                      phi = (phi + 2 * Math.PI) % (2 * Math.PI);
                      theta = 2 * Math.PI * (onMouseDownPosition.y - event.clientY) / 400 + onMouseDownTheta;
                      var epsilon = 1e-12; // Prevents spinnging from getting stuck
                      theta = Math.max(Math.min(Math.PI - epsilon, theta), epsilon);
                
                      update_camera_position();
                    }
                    render();
                  } else {
                    container.style.cursor = "pointer";
                  }
                }
                
                function onDocumentMouseUp(event) {
                  event.preventDefault();
                
                  isMouseDown = false;
                  container.style.cursor = "pointer";
                
                  if (autoRescale) {
                      ScaleInView();
                      render();
                  }
                  positionAxes();
                  render();
                  positionticknums();
                }
                
                // Bind Mouse events
                container.addEventListener('mousemove', onDocumentMouseMove, false);
                container.addEventListener('mousedown', onDocumentMouseDown, false);
                container.addEventListener('mouseup', onDocumentMouseUp, false);
                onMouseDownPosition = new THREE.Vector2();
                var autoRescale = true;
                
                update_camera_position();
                positionAxes();
                render(); // Rendering twice updates camera.matrixWorldInverse so that ScaleInView works properly
                ScaleInView();
                render();     
                positionticknums();	
				

				
			break;
			
			case 'List':
				var copy = Object.assign({}, params);
				
				var mess = [];
				
				func.args.forEach(function(el) {
					mess.push(interpretate(el, func, copy, mesh));
				});
				
				return mess;
			
			break;
	
			case 'Style':
				var copy = Object.assign({}, params);
				
				func.args.forEach(function(el) {
					interpretate(el, func, copy, mesh);
				});
			
			break;		
	
			case 'Annotation':
				
				func.args.forEach(function(el) {
					interpretate(el, func, params, mesh);
				});
			
			break;		
			
			case 'GraphicsGroup':
				var group = new THREE.Group();
				var copy = Object.assign({}, params);
				
				func.args.forEach(function(el) {
					interpretate(el, func, copy, group);
				});
				
				mesh.add(group);
			
			break;
			
			case 'RGBColor':
				if (func.args.length !== 3) console.error( "RGB values should be triple!");
			
				var r = Math.round(255*interpretate(func.args[0]));
				var g = Math.round(255*interpretate(func.args[1]));
				var b = Math.round(255*interpretate(func.args[2]));
				
				params.color = new THREE.Color("rgb("+r+","+g+","+b+")");			
			
			break;
			
			case 'Opacity':
				var o = interpretate(func.args[0]);
				if (typeof o !== 'number') console.error( "Opacity must have number value!");
				console.log(o);
				params.opacity = o;
			
			break;
			
			case 'Thickness':
				//params.thickness = 
			break;
			
			case 'Arrowheads':
			
			break;
			
			case 'Arrow':
				interpretate(func.args[0], func, params, mesh);
				
			break;
			
			case 'Tube':
				var arr = interpretate(func.args[0]);
				if (arr.length ==  1) arr = arr[0];
				if (arr.length !== 2) console.error( "Tube must have 2 vectors!");
				
				var points = [new THREE.Vector4(...arr[0], 1),
							new THREE.Vector4(...arr[1], 1)];
				
				points.forEach(function(p) {
					p = p.applyMatrix4(params.matrix);
				});
	
				var origin = points[0].clone();
				var dir = points[1].add(points[0].negate());
				
				var arrowHelper = new THREE.ArrowHelper(dir.normalize(), origin, dir.length(), params.color);
				mesh.add(arrowHelper);
				arrowHelper.line.material.linewidth = params.thickness;
			break;
	
			case 'Sphere':
				var radius = 1;
				if (func.args.length > 1) radius = func.args[1];
				
				var material = new THREE.MeshLambertMaterial({
					color:params.color,
					transparent:false,
					opacity:params.opacity,
				});				
				
				function addSphere(cr) {
					var origin = new THREE.Vector4(...cr, 1);
					var geometry = new THREE.SphereGeometry( radius, 20, 20 );
					var sphere = new THREE.Mesh( geometry, material );
					
					sphere.position.x = origin.x; 
					sphere.position.y = origin.y; 
					sphere.position.z = origin.z;
		
					mesh.add( sphere );
					geometry.dispose(); 			
				}	
	
				var list = interpretate(func.args[0]);
				
				if (list.length == 1) list=list[0];
				if (list.length == 1) list=list[0];
				
				if (list.length == 3) {
					addSphere(list);
				} else if (list.length > 3) {
					list.forEach(function(el) {
						addSphere(el);		
					});	
				} else {
					console.log(list);
					console.error( "List of coords. for sphere object is less 1");
				}
				
				material.dispose();
				
				break;
				
			case 'Cuboid':
				//if (params.hasOwnProperty('geometry')) {
				//	var points = [new THREE.Vector4(...interpretate(func.args[0]), 1),
				//				new THREE.Vector4(...interpretate(func.args[1]), 1)];				
				//}
				console.log("Cuboid");
				var points, diff, origin;
				
				if (func.args.length == 2) {
				
					points = [new THREE.Vector4(...interpretate(func.args[0]), 1),
							new THREE.Vector4(...interpretate(func.args[1]), 1)];			
				
					origin = points[0].clone().add(points[1].clone()).divideScalar(2);
					diff = points[0].clone().add(points[1].clone().negate());
				
				} else if (func.args.length == 1) {
					p = interpretate(func.args[0]);
					origin = new THREE.Vector4(...p, 1);
					diff = new THREE.Vector4(1,1,1, 1);
					
					//shift it
					origin.add(diff.clone().divideScalar(2))
				} else {
					console.error( "Expected 2 or 1 arguments");
				}
				
				
				var geometry = new THREE.BoxGeometry(diff.x, diff.y, diff.z);
				var material = new THREE.MeshLambertMaterial({
					color:params.color,
					transparent:true,
					opacity:params.opacity,
					depthWrite: true
				});
				//material.side = THREE.DoubleSide;
				
				var cube = new THREE.Mesh( geometry, material );
				
				//var tr = new THREE.Matrix4();
				//	tr.makeTranslation(origin.x,origin.y,origin.z);
				
				//cube.applyMatrix(params.matrix.clone().multiply(tr));
				
				
				
				cube.position.x = origin.x; 
				cube.position.y = origin.y; 
				cube.position.z = origin.z; 
				
				
				mesh.add(cube);
				
				geometry.dispose();
				material.dispose();
			
			break;
			
			case 'Rule':
				switch(func.args[0]) {
					
					
				}
			break;
			
			case 'Center':
				return 'Center';
			break;
			
			case 'Tetrahedron':
				var points = interpretate(func.args[0]);
				var faces = [];
				console.log("Points of tetrahedron:");
				console.log(points);
				
				faces.push([points[0], points[1], points[2]]);
				faces.push([points[0], points[1], points[3]]);
				faces.push([points[1], points[2], points[3]]);
				faces.push([points[0], points[3], points[2]]);

				var fake = ["List"];
				
				var listVert = function(cord) {
						return([
								"List",
								cord[0],
								cord[1],
								cord[2]
								]);
				}
					
				faces.forEach(function(fs) {

					var struc = [
									"Polygon",
									[
										"List",
										listVert(fs[0]),
										listVert(fs[1]),
										listVert(fs[2])
									]
								];
					fake.push(struc);
				});
				console.log(fake);
				interpretate(fake, func, params, mesh);
			break;
					
			
			case 'GeometricTransformation':
				
				var group = new THREE.Group();
				//Если center, то наверное надо приметь matrix 
				//к каждому объекту относительно родительской группы.
				var p = interpretate(func.args[1]);
				var centering = false;
				var centrans = [];
				
				
				
				
				if (p.length === 1) { p = p[0]; }
				if (p.length === 1) { p = p[0]; }
				else if (p.length === 2) {
					console.log(p);
					if (p[1] === 'Center') {
						
						centering = true;
					} else {
						console.log("NON CENTERING ISSUE!!!");
						console.log(p);
						centrans = p[1];
						console.log("???");
					}
					//return;
					p = p[0];
				}
				
				if (p.length === 3) {
					
				
					if (typeof p[0] === 'number') {
						var dir = p;
						var matrix = new THREE.Matrix4().makeTranslation(...dir,1);
					} else {
						
						//make it like Matrix4
						p.forEach(function(el) {
							el.push(0);
						});
						p.push([0,0,0,1]);
						
			
						var matrix = new THREE.Matrix4();
						console.log("Apply matrix to group::");
						matrix.set(...flatten(p));
					}
				} else {
					console.log(p);
					console.error( "Unexpected length matrix: :: " + p);
				}
				
				
				//Backup of params
				var copy = Object.assign({}, params);	
				interpretate(func.args[0], func, copy, group);
				
				console.log(matrix);
				
				
				if (centering || centrans.length > 0) {
					console.log("::CENTER::");
					var bbox = new THREE.Box3().setFromObject(group);
					console.log(bbox);
					var center = new THREE.Vector3().addVectors(bbox.max,bbox.min).divideScalar(2);
					if (centrans.length > 0) {
						console.log("CENTRANS");
						center = center.fromArray(centrans);
					}
					console.log(center);
					
					var	translate = new THREE.Matrix4().makeTranslation(-center.x,-center.y,-center.z,1);
					group.applyMatrix(translate);
					group.applyMatrix(matrix);
						translate = new THREE.Matrix4().makeTranslation(center.x,center.y,center.z,1);
					group.applyMatrix(translate);
				} else {
					group.applyMatrix(matrix);
				}
				
				mesh.add(group);
				
			break;
	
			case 'GraphicsComplex':			
				var copy = Object.assign({}, params);
				
				copy.geometry = new THREE.Geometry();
				
				
				interpretate(func.args[0]).forEach(function(el) {
					if (typeof el[0] !== 'number') console.error( "not a triple of number"+el);
					copy.geometry.vertices.push(
						new THREE.Vector3(el[0],el[1],el[2])
					);
				});
	
				
				var group = new THREE.Group();
				
				
				interpretate(func.args[1], func, copy, mesh);
				
				mesh.add(group);	
				copy.geometry.dispose();
				
			break;
	
			case 'Polygon':
				if (params.hasOwnProperty('geometry')) {
					var geometry = params.geometry.clone();
	
	
					var createFace = function(c) {
						
						switch(c.length) {
							case 3:
								geometry.faces.push(new THREE.Face3(c[0]-1,c[1]-1,c[2]-1));
							break;
							
							case 4:
								geometry.faces.push(new THREE.Face3(c[0]-1,c[1]-1,c[2]-1));
								geometry.faces.push(new THREE.Face3(c[0]-1,c[2]-1,c[3]-1));
							break;
							
							case 5:
								geometry.faces.push(new THREE.Face3(c[0]-1,c[1]-1,c[4]-1));
								geometry.faces.push(new THREE.Face3(c[1]-1,c[2]-1,c[3]-1));
								geometry.faces.push(new THREE.Face3(c[1]-1,c[3]-1,c[4]-1));
							break;
							
							default:
								console.log(c);
								console.log(c.length);
								console.error( "Cant produce complex polygons! at"+c);
							
						}
					}
					
					var a = interpretate(func.args[0]);
					if (a.length === 1) {
						a = a[0];
					}
					
					if (typeof a[0] === 'number') {
						console.log("Create single face");
						createFace(a);
					} else {
						console.log("Create multiple face");
						console.log(a);
						a.forEach(function(el) {
						
							createFace(el);
						});
					}
				} else {
					
					var geometry = new THREE.Geometry();
					var points = interpretate(func.args[0]);
				
					points.forEach(function(el) {
						if (typeof el[0] !== 'number') console.error( "not a triple of number"+el);
						geometry.vertices.push(
							new THREE.Vector3(el[0],el[1],el[2])
						);					
					});
					
					console.log("points");
					console.log(points);
					
					switch(points.length) {
						case 3:
							geometry.faces.push(new THREE.Face3(0,1,2));
						break;
						
						case 4:
							geometry.faces.push(new THREE.Face3(0,1,2));
							geometry.faces.push(new THREE.Face3(0,2,3));
						break;
						
						case 5:
							geometry.faces.push(new THREE.Face3(0,1,4));
							geometry.faces.push(new THREE.Face3(1,2,3));
							geometry.faces.push(new THREE.Face3(1,3,4));
						break;
						
						default:
							console.log(points);
							console.error( "Cant build complex polygon ::");
					}
					
				}
				
				var material = new THREE.MeshLambertMaterial({
					color:params.color,
					transparent: params.opacity < 0.9? true : false,
					opacity:params.opacity,

					//depthTest: false
					//depthWrite: false
				});
				console.log(params.opacity);
				material.side = THREE.DoubleSide;
				
				geometry.computeFaceNormals();
				//complex.computeVertexNormals();
				var poly = new THREE.Mesh(geometry, material);
				
				//poly.frustumCulled = false;
				mesh.add(poly);
				material.dispose();
				
			break;
			
			case 'GrayLevel':
			
			break;
			
			case 'EdgeForm':
			
			break;	
	
			case 'Specularity':
				
			break;
			
			case 'Text':
			
			break;
			
			case 'Line':
				if (params.hasOwnProperty('geometry')) {
					var geometry = new THREE.Geometry();
					
					var points = interpretate(func.args[0]);
					points.forEach(function(el) {
						geometry.vertices.push(new THREE.Vector3().copy(params.geometry.vertices[el-1]));
					});
					
					var material = new THREE.LineBasicMaterial( { linewidth: params.thickness, color: params.edgecolor } );
					var line = new THREE.Line( geometry, material );
					
					line.material.polygonOffset = true;
					line.material.polygonOffsetFactor = 1;
					line.material.polygonOffsetUnits = 1;
					
					mesh.add(line);
					
					geometry.dispose();
					material.dispose();
				} else {
					var arr = interpretate(func.args[0]);
					if (arr.length ==  1) arr = arr[0];
					//if (arr.length !== 2) console.error( "Tube must have 2 vectors!");
					console.log("points: "+arr.length);
				
					var points = [];
					arr.forEach(function(p) {
						points.push(new THREE.Vector4(...p, 1));
					});
					//new THREE.Vector4(...arr[0], 1)
					
					points.forEach(function(p) {
						p = p.applyMatrix4(params.matrix);
					});
		
					const geometry = new THREE.BufferGeometry().setFromPoints( points );
					const material = new THREE.LineBasicMaterial({
						color: params.edgecolor,
						linewidth: params.thickness
					});
					
					mesh.add(new THREE.Line( geometry, material ));
					
				}
			break;
			
			default:
				console.error( "Undefined function : "+func.name);
				
			break;
			
		}
		
		return(undefined);
		
	}


