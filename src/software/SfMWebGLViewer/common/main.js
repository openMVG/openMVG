// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The canvas containing the webgl object 
var divCanvas; // The div containing the canvas object 
var canvas;    // The canvas itself 

// The render context (GL context, shaders)
var renderContext ; 

// Some state variables of the interface 
var optionEnableCameraBtn ;
var cameraEnabled ;
var optionEnablePointCloudBtn ;
var pcloudEnabled ; 

// The trackball object 
var trackball ; 
var optionEnableTrackballBtn ; 

// The camera 
var camera ; 

// The pointcloud
var pointcloud ; 

function init()
{
  divCanvas = document.getElementById("divViewer");
  canvas = document.getElementById("viewerCanvas");

  optionEnableCameraBtn = document.getElementById("divShowCam");
  optionEnablePointCloudBtn = document.getElementById("divShowPCloud");
  optionEnableTrackballBtn = document.getElementById("divShowTrackball");

  cameraEnabled = true ; 
  pcloudEnabled = true ; 

  renderContext = new RenderContext( canvas ) ; 
  trackball = new Trackball( 128 , renderContext ) ; 

  var camPos = Vector.create( 0 , 0 , -10 ) ;
  var camDst = Vector.create( 0.0 , 0.0 , 0.0 ) ; 
  var camUp = Vector.create( 0.0 , -1.0 , 0.0 ) ; 

  camera = new PerspectiveCamera( camPos , camDst , camUp , 75.0 , 0.1 , 1000.0 , 640 , 480 ) ;  

  /* Center scene on the pcloud */
  var bs = ComputeBoundingSphere( modelPos ) ;
  camera.FitBoundingSphere( bs , 1.2 * bs[3] ) ; 
  trackball.setPosition( bs ) ;
  trackball.setRadius( bs[3] ) ; 

  setEventListeners();

  pointcloud = new PointCloud( modelPos , undefined , modelCol , renderContext ) ; 

  // Set scene for the first time 
  update() ; 
  resizeGLCanvas();
}


function toggleEnableCameras()
{
  if( cameraEnabled )
  {
    optionEnableCameraBtn.className = optionEnableCameraBtn.className.replace( 'active' , '' ) ;
    optionEnableCameraBtn.classList ? optionEnableCameraBtn.classList.add('inactive') : optionEnableCameraBtn.className += ' inactive';
    
    cameraEnabled = false ; 
  }
  else 
  {
    optionEnableCameraBtn.className = optionEnableCameraBtn.className.replace( 'inactive' , '' ) ;
    optionEnableCameraBtn.classList ? optionEnableCameraBtn.classList.add('active') : optionEnableCameraBtn.className += ' active';

    cameraEnabled = true ; 
  }
  update();
}

function toggleEnablePointCloud()
{
  if( pcloudEnabled )
  {
    optionEnablePointCloudBtn.className = optionEnablePointCloudBtn.className.replace( 'active' , '' ) ;
    optionEnablePointCloudBtn.classList ? optionEnablePointCloudBtn.classList.add('inactive') : optionEnablePointCloudBtn.className += ' inactive';
    
    pcloudEnabled = false ; 
  }
  else 
  {
    optionEnablePointCloudBtn.className = optionEnablePointCloudBtn.className.replace( 'inactive' , '' ) ;
    optionEnablePointCloudBtn.classList ? optionEnablePointCloudBtn.classList.add('active') : optionEnablePointCloudBtn.className += ' active';

    pcloudEnabled = true ; 
  }
  update(); 
}

function toggleEnableTrackball()
{
  if( trackball.isVisible() )
  {
    RemoveClass( optionEnableTrackballBtn , 'active' ) ;
    AddClass( optionEnableTrackballBtn , 'inactive' ) ; 

    trackball.setVisible( false ) ;
  }
  else
  {
    RemoveClass( optionEnableTrackballBtn , 'inactive' ) ;
    AddClass( optionEnableTrackballBtn , 'active' ) ; 

    trackball.setVisible( true ) ; 
  }
  update(); 
}

function setEventListeners()
{
  // Buttons 
  optionEnableCameraBtn.addEventListener( "click" , toggleEnableCameras ) ;
  optionEnablePointCloudBtn.addEventListener( "click" , toggleEnablePointCloud ) ; 
  optionEnableTrackballBtn.addEventListener( "click" , toggleEnableTrackball ) ; 

  // Window resizing 
  window.addEventListener( 'resize' , resizeGLCanvas ) ; 
}

function resizeGLCanvas()
{
  var w = divCanvas.clientWidth ;
  var h = divCanvas.clientHeight ;

  if( Math.abs( canvas.width - w ) > 10 || Math.abs( canvas.height - h ) > 10 )
  {
    canvas.width = w ;
    canvas.height = h ;

    camera.m_width = w ;
    camera.m_height = h ; 
  }
  draw();
}

function update()
{
  // Update internal matrices 
  renderContext.setCurrentViewMatrix( camera.GetLookAtMatrix() );
  renderContext.setCurrentProjectionMatrix( camera.GetProjectionMatrix() ) ;

  // Finally draw the scene 
  draw() ; 
}

function draw()
{
  var gl = renderContext.getGLContext() ; 

  gl.viewport( 0 , 0 , canvas.width , canvas.height ) ; 
  gl.clear( gl.COLOR_BUFFER_BIT , gl.DEPTH_BUFFER_BIT ) ;

  // Draw the trackball   
  trackball.draw( renderContext ) ; 

  // Draw the scene 
  pointcloud.draw( renderContext ) ; 

  // Draw the cameras 
}
