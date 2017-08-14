// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The canvas containing the webgl object 
var divCanvas; // The div containing the canvas object 
var canvas;    // The canvas itself 

// The render context (GL context, shaders)
var renderContext; 

// Some state variables of the interface 
var optionEnableCameraBtn;
var cameraEnabled;

// The trackball object 
var trackball; 
var optionEnableTrackballBtn; 

// The camera 
var camera; 

// The pointcloud
var pointcloud; 
var optionEnablePointCloudBtn;
var pointSizeSlider; 
var optionViewDependentPointSizeBtn; 
var optionBackFaceCullingBtn; 

// The gizmos 
var cameraGizmos;  

var optionResetBtn;
var cameraScaleSlider; 

var mouseEnterPosition = {};
var mouseLastPosition = {}; 
var mouseIsClicked = false;

// Main entry point 
function init()
{
  divCanvas = document.getElementById("divViewer");
  canvas = document.getElementById("viewerCanvas");

  // Btn
  optionEnableCameraBtn = document.getElementById("divShowCam");
  optionEnablePointCloudBtn = document.getElementById("divShowPCloud");
  optionEnableTrackballBtn = document.getElementById("divShowTrackball");
  optionResetBtn = document.getElementById( "divResetButton" );
  optionBackFaceCullingBtn = document.getElementById( "divBackFaceCullingButton" ); 
  optionViewDependentPointSizeBtn = document.getElementById( "divViewDependendPointSizeBtn" ); 

  // Sliders 
  cameraScaleSlider = document.getElementById( "cameraScaleSlider" );
  pointSizeSlider = document.getElementById( "pointSizeSlider" ); 
  
  cameraEnabled = true; 
  pcloudEnabled = true; 

  renderContext = new RenderContext( canvas ); 
  trackball = new Trackball( 128 , renderContext ); 

  var camPos = Vector.create( 0 , 0 , -10 );
  var camDst = Vector.create( 0.0 , 0.0 , 0.0 ); 
  var camUp = Vector.create( 0.0 , -1.0 , 0.0 ); 

  camera = new PerspectiveCamera( camPos , camDst , camUp , 75.0 , 0.01 , 1000.0 , 640 , 480 );  

  cameraGizmos = new Array( cameraPos.length );
  for (var i = 0; i < cameraPos.length; ++i )
  {
    cameraGizmos[ i ] = new CameraGizmo( cameraPos[i] , cameraImagePlanes[i] , undefined , renderContext );
  }

  /* Center scene on the pcloud */
  var bs = ComputeBoundingSphere( modelPos );
  camera.FitBoundingSphere( bs , bs[3] ); 
  trackball.setPosition( bs );
  var d = Vector.norm( Vector.sub( camera.m_pos , camera.m_dir ) );
  var arad = DegToRad( camera.m_fov ); 
  var new_rad = d * Math.sin( arad * 0.5 ); 
  trackball.setRadius( bs[3] * 0.8 ); 

  setEventListeners();

  pointcloud = new PointCloud( modelPos , modelNor , modelCol , renderContext ); 

  // Set scene for the first time 
  update(); 
  resizeGLCanvas();
}

// Handle visibility on the camera 
function toggleEnableCameras()
{
  if (cameraEnabled )
  {
    optionEnableCameraBtn.className = optionEnableCameraBtn.className.replace( 'active' , '' );
    optionEnableCameraBtn.classList ? optionEnableCameraBtn.classList.add('inactive') : optionEnableCameraBtn.className += ' inactive';
    
    for (var i = 0; i < cameraGizmos.length; ++i )
    {
      cameraGizmos[ i ].setVisible( false ); 
    }
    cameraEnabled = false; 
  }
  else 
  {
    optionEnableCameraBtn.className = optionEnableCameraBtn.className.replace( 'inactive' , '' );
    optionEnableCameraBtn.classList ? optionEnableCameraBtn.classList.add('active') : optionEnableCameraBtn.className += ' active';

    for (var i = 0; i < cameraGizmos.length; ++i )
    {
      cameraGizmos[ i ].setVisible( true ); 
    }
    cameraEnabled = true; 
  }
  update();
}

// Handle visibility on the point cloud 
function toggleEnablePointCloud()
{
  if (pointcloud.isVisible() )
  {
    optionEnablePointCloudBtn.className = optionEnablePointCloudBtn.className.replace( 'active' , '' );
    optionEnablePointCloudBtn.classList ? optionEnablePointCloudBtn.classList.add('inactive') : optionEnablePointCloudBtn.className += ' inactive';

    pointcloud.setVisible( false );
  }
  else 
  {
    optionEnablePointCloudBtn.className = optionEnablePointCloudBtn.className.replace( 'inactive' , '' );
    optionEnablePointCloudBtn.classList ? optionEnablePointCloudBtn.classList.add('active') : optionEnablePointCloudBtn.className += ' active';

    pointcloud.setVisible( true ); 
  }
  update(); 
}

// Handle visibility on the trackball
function toggleEnableTrackball()
{
  if (trackball.isVisible() )
  {
    RemoveClass( optionEnableTrackballBtn , 'active' );
    AddClass( optionEnableTrackballBtn , 'inactive' ); 

    trackball.setVisible( false );
  }
  else
  {
    RemoveClass( optionEnableTrackballBtn , 'inactive' );
    AddClass( optionEnableTrackballBtn , 'active' ); 

    trackball.setVisible( true ); 
  }
  update(); 
}

function toggleBackfaceCulling()
{
  if (renderContext.isBackfaceCullingActive() )
  {
    RemoveClass( optionBackFaceCullingBtn , 'active' );
    AddClass( optionBackFaceCullingBtn , 'inactive' );
    
    renderContext.setBackFaceCulling( false ); 
  }
  else 
  {
    RemoveClass( optionBackFaceCullingBtn , 'inactive' );
    AddClass( optionBackFaceCullingBtn , 'active' ); 

    renderContext.setBackFaceCulling( true ); 
  }
  update(); 
}

function toggleViewDependentPointSize()
{
  if (renderContext.isViewDependentPointSizeActive() )
  {
    RemoveClass( optionViewDependentPointSizeBtn , 'active' );
    AddClass( optionViewDependentPointSizeBtn , 'inactive' ); 
    
    renderContext.setViewDependentPointSize( false ); 
  }
  else 
  {
    RemoveClass( optionViewDependentPointSizeBtn , 'inactive' );
    AddClass( optionViewDependentPointSizeBtn , 'active' );

    renderContext.setViewDependentPointSize( true ); 
  }
  update(); 
}

// Handle mouse click 
function onMouseDown( e )
{
  mouseEnterPosition.x = e.clientX - canvas.offsetLeft;
  mouseEnterPosition.y = e.clientY - canvas.offsetTop;

  mouseLastPosition.x = e.clientX - canvas.offsetLeft;
  mouseLastPosition.y = e.clientY - canvas.offsetTop;

  if (e.which === 1 )
  {
    mouseIsClicked = true; 
  } 
}

// Handle mouse unclick 
function onMouseUp( e ) 
{
  if (e.which === 1 )
  {
    mouseIsClicked = false; 
  }
}

// Handle mouse drag 
function onMouseMove( e )
{
  if (e.button == 0 && mouseIsClicked ) // Left click 
  {
    var mouseCurX = e.clientX - canvas.offsetLeft;
    var mouseCurY = e.clientY - canvas.offsetTop;

    if (e.shiftKey )
    { 
      var dx = mouseCurX - mouseLastPosition.x;
      var dy = mouseCurY - mouseLastPosition.y;

      // Zoom camera 
      camera.zoom( dy / 20.0 );   

      // Update trackball 
      var d = Vector.norm( Vector.sub( camera.m_pos , camera.m_dir ) );
      var arad = DegToRad( camera.m_fov ); 
      var new_rad = d * Math.sin( arad * 0.5 ); 
      trackball.setRadius( new_rad * 0.8 ); 
    }
    else if (e.altKey )
    {
      // Pan 
      var p_old = camera.pointOnPlane( mouseLastPosition.x , canvas.height - mouseLastPosition.y );
      var p_new = camera.pointOnPlane( mouseCurX , canvas.height - mouseCurY );

      var d = Vector.sub( p_new , p_old ); 

      pointcloud.translate( d ); 
      for (var i = 0; i < cameraGizmos.length; ++i )
      {
        cameraGizmos[ i ].translate( d ); 
      }  
    }
    else 
    {
      var p_old = camera.pointOnSphere( mouseLastPosition.x , canvas.height - mouseLastPosition.y , trackball.getRadius() );
      var p_new = camera.pointOnSphere( mouseCurX , canvas.height - mouseCurY , trackball.getRadius() );

      // Axis of rotation 
      var d1 = Vector.sub( p_new , camera.m_dir );
      var d2 = Vector.sub( p_old , camera.m_dir ); 
      var axis = Vector.cross( d1 , d2 );
      
      d1 = Vector.normalize( d1 );
      d2 = Vector.normalize( d2 ); 
      var d = Vector.norm( Vector.sub( p_old , p_new ) );
      // Angle of rotation 
      var angle = 0.75 * d / trackball.getRadius(); 

      if (angle > 0.0 )
      {
        var q = new Quaternion( 0.0 , 0.0 , 0.0 , 0.0 );
        q.setFromAxisAngle( axis , - angle ); 


        var neg_dir = new Vector.negate( camera.m_dir ); 
        var dqv = new DualQuaternion();
        dqv.setFromTranslationVector( camera.m_dir );

        var invdqv = new DualQuaternion();
        invdqv.setFromTranslationVector( neg_dir ); 
        
        var dqq = new DualQuaternion(); 
        dqq.setFromRotationQuaternion( q ); 

        // Update trackball
        trackball.m_orient = dqv.mul( dqq.mul( invdqv.mul( trackball.m_orient ) ) );

        // Update point cloud 
        pointcloud.rotate( q , camera.m_dir ); 
        for (var i = 0; i < cameraGizmos.length; ++i )
        {
          cameraGizmos[ i ].rotate( q , camera.m_dir ); 
        }
      }
      
      // Rotate 
    }

    mouseLastPosition.x = mouseCurX;
    mouseLastPosition.y = mouseCurY;

    update( );
  }
}

function onSliderChange()
{
  var scaleValue = cameraScaleSlider.value; 
  for (var i = 0; i < cameraGizmos.length; ++i )
  {
    cameraGizmos[i].setScale( scaleValue ); 
  }
  update(); 
}

function onPointSizeSliderChange()
{
  var pointSizeValue = pointSizeSlider.value; 
  pointcloud.setPointSize( pointSizeValue ); 
  update(); 
}

// Setup event listener 
function setEventListeners()
{
  // Buttons 
  optionEnableCameraBtn.addEventListener( "click" , toggleEnableCameras );
  optionEnablePointCloudBtn.addEventListener( "click" , toggleEnablePointCloud ); 
  optionEnableTrackballBtn.addEventListener( "click" , toggleEnableTrackball ); 
  optionResetBtn.addEventListener( "click" , resetView ); 
  optionBackFaceCullingBtn.addEventListener( "click" , toggleBackfaceCulling ); 
  divViewDependendPointSizeBtn.addEventListener( "click" , toggleViewDependentPointSize ); 

  // Canvas GL Window
  canvas.addEventListener( "mousedown" , onMouseDown );
  canvas.addEventListener( "mousemove" , onMouseMove );
  canvas.addEventListener( "mouseup" , onMouseUp );

  // Sliders 
  cameraScaleSlider.addEventListener( "change" , onSliderChange ); 
  cameraScaleSlider.addEventListener( "input" , onSliderChange ); 
  pointSizeSlider.addEventListener( "change" , onPointSizeSliderChange );
  pointSizeSlider.addEventListener( "input" , onPointSizeSliderChange ); 
  

  // Window resizing 
  window.addEventListener( 'resize' , resizeGLCanvas ); 
}

// Function called when user resize the window 
function resizeGLCanvas()
{
  var w = divCanvas.clientWidth;
  var h = divCanvas.clientHeight;

  if (Math.abs( canvas.width - w ) > 10 || Math.abs( canvas.height - h ) > 10 )
  {
    canvas.width = w;
    canvas.height = h;

    camera.m_width = w;
    camera.m_height = h; 
  }
  update();
}

/* Update all rendering infos before rendering (matrices) */
function update()
{
  // Update internal matrices 
  renderContext.setCurrentViewMatrix( camera.GetLookAtMatrix() );
  renderContext.setCurrentProjectionMatrix( camera.GetProjectionMatrix() );
  
  // Set current camera 
  renderContext.setCurrentCamera( camera ); 

  // Finally draw the scene 
  draw(); 
}

/* Draw the scene */
function draw()
{
  var gl = renderContext.getGLContext(); 

  gl.viewport( 0 , 0 , canvas.width , canvas.height ); 
  gl.clear( gl.COLOR_BUFFER_BIT , gl.DEPTH_BUFFER_BIT );

  // Draw the trackball   
  trackball.draw( renderContext ); 

  // Draw the scene 
  pointcloud.draw( renderContext ); 

  // Draw the cameras 
  for (var i = 0; i < cameraGizmos.length; ++i )
  {
    cameraGizmos[ i ].draw( renderContext ); 
  }
}

/* Reset the view on the default orientation */
function resetView()
{
  trackball.reset();
  pointcloud.reset();
  for (var i = 0; i < cameraGizmos.length; ++i )
  {
    cameraGizmos[ i ].reset(); 
  }

  var camPos = Vector.create( 0 , 0 , -10 );
  var camDst = Vector.create( 0.0 , 0.0 , 0.0 ); 
  var camUp = Vector.create( 0.0 , -1.0 , 0.0 ); 

  camera = new PerspectiveCamera( camPos , camDst , camUp , 75.0 , 0.01 , 1000.0 , canvas.width , canvas.height );  

  /* Center scene on the pcloud */
  var bs = ComputeBoundingSphere( modelPos );
  camera.FitBoundingSphere( bs , bs[3] ); 
  trackball.setPosition( bs );
  var d = Vector.norm( Vector.sub( camera.m_pos , camera.m_dir ) );
  var arad = DegToRad( camera.m_fov ); 
  var new_rad = d * Math.sin( arad * 0.5 ); 
  trackball.setRadius( bs[3] * 0.8 ); 
  
  // Default values for interface
  pointSizeSlider.value = 1.0; 
  cameraScaleSlider.value = 1.0; 

  RemoveClass( optionBackFaceCullingBtn , 'inactive' );
  RemoveClass( optionBackFaceCullingBtn , 'active' );
  AddClass( optionBackFaceCullingBtn , 'active' ); 
  renderContext.setBackFaceCulling( true ); 

  RemoveClass( optionViewDependentPointSizeBtn , 'inactive' );
  RemoveClass( optionViewDependentPointSizeBtn , 'active' );
  AddClass( optionViewDependentPointSizeBtn , 'active' ); 
  renderContext.setViewDependentPointSize( true ); 

  update();
}