var io = require('socket.io-client');
var ss = require('socket.io-stream');
const fs = require('fs');

var socket = io();
var stream = ss.createStream({
  highWaterMark: 1024,
  objectMode: true
});

ss(socket).emit('setup-stream', stream);

main();

function main() {

  const canvas = document.querySelector('#glcanvas');
  const gl = canvas.getContext('webgl');

  if (!gl) {
    alert('Unable to initialize WebGL. Your browser or machine may not support it.');
    return;
  }

  const vsSource = `
    attribute vec4 aVertexPosition;

    uniform mat4 uModelViewMatrix;
    uniform mat4 uProjectionMatrix;

    void main() {
      gl_Position = uProjectionMatrix * uModelViewMatrix * aVertexPosition;
      gl_PointSize = 1.0;
    }
  `;

  const fsSource = `
    void main() {
      gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
    }
  `;

  const shaderProgram = initShaderProgram(gl, vsSource, fsSource);

  const programInfo = {
    program: shaderProgram,
    attribLocations: {
      vertexPosition: gl.getAttribLocation(shaderProgram, 'aVertexPosition'),
    },
    uniformLocations: {
      projectionMatrix: gl.getUniformLocation(shaderProgram, 'uProjectionMatrix'),
      modelViewMatrix: gl.getUniformLocation(shaderProgram, 'uModelViewMatrix'),
    },
  };

  const vsSource1 = `
    attribute vec4 aVertexPosition;
    attribute vec2 aTexCoord;
    varying vec2 vTexCoord;

    void main() {
      gl_Position = aVertexPosition;
      vTexCoord = aTexCoord;
    }
  `;

  const fsSource1 = `
    precision mediump float;
    varying vec2 vTexCoord;
    uniform sampler2D uTexture;

    void main() {
      gl_FragColor = texture2D(uTexture, vTexCoord);
    }
  `;

  const shaderProgram1 = initShaderProgram(gl, vsSource1, fsSource1);

  const programInfo1 = {
    program: shaderProgram1,
    attribLocations: {
      vertexPosition: gl.getAttribLocation(shaderProgram1, 'aVertexPosition'),
      textureCoord: gl.getAttribLocation(shaderProgram1, 'aTexCoord'),
    },
    uniformLocations: {
      textureSampler: gl.getUniformLocation(shaderProgram1, 'uTexture')
    },
  };

  const buffers = initBuffers(gl);

  stream.on('readable', function() {

    var positions = new Float32Array( 2 * 150 );
    var tmp = stream.read(1);//maybe we can get a Float32Array right here.
    for (var i = 0; i < 2 * 150; i++) {
      positions[i] = tmp[i];
    }

    gl.bindBuffer(gl.ARRAY_BUFFER, buffers.position);
    gl.bufferData(gl.ARRAY_BUFFER,
                  positions,
                  gl.DYNAMIC_DRAW);

    drawScene(gl, programInfo, programInfo1, buffers);
  });

}

function initBuffers(gl) {

  const positionBuffer = gl.createBuffer();
  const quadvbo= gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, quadvbo);
  const quad = [
     1.0,  1.0, 1.0, 1.0,
    -1.0,  1.0, 0.0, 1.0,
     1.0, -1.0, 1.0, 0.0,
    -1.0, -1.0, 0.0, 0.0
  ];
  gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(quad),gl.STATIC_DRAW);
  const frameBuffer = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuffer);

  texture = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, texture);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA,
                gl.canvas.clientWidth, gl.canvas.clientHeight, 0,
                gl.RGBA, gl.UNSIGNED_BYTE, null);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);

  /*const rbo = gl.createRenderbuffer();
  gl.bindRenderbuffer(gl.RENDERBUFFER, rbo);
  gl.renderbufferStorage(gl.RENDERBUFFER, gl.RGBA4, gl.canvas.clientWidth, gl.canvas.clientHeight);
  gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.RENDERBUFFER, rbo);*/
  if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) {
   alert("Framebuffer setup failed.");
  }

  return {
    position: positionBuffer,
    fb: frameBuffer,
    quadvbo: quadvbo,
    texture: texture
  };
}

function drawScene(gl, programInfo, programInfo1, buffers) {

  gl.bindFramebuffer(gl.FRAMEBUFFER, buffers.fb);
  //gl.clearDepth(1.0);
  //gl.enable(gl.DEPTH_TEST);
  //gl.depthFunc(gl.LEQUAL);
  //gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

  const fieldOfView = 60 * Math.PI / 180;
  const aspect = gl.canvas.clientWidth / gl.canvas.clientHeight;
  const zNear = 0.1;
  const zFar = 400.0;

  const projectionMatrix = mat4.create();
  mat4.perspective(projectionMatrix,
                   fieldOfView,
                   aspect,
                   zNear,
                   zFar);

  const modelViewMatrix = mat4.create();
  mat4.translate(modelViewMatrix,
                 modelViewMatrix,
                 [-50.0, -50.0, -250.0]);

  gl.bindBuffer(gl.ARRAY_BUFFER, buffers.position);
  gl.vertexAttribPointer(
        programInfo.attribLocations.vertexPosition,
        2,
        gl.FLOAT,
        false,
        0,
        0);
  gl.enableVertexAttribArray(programInfo.attribLocations.vertexPosition);

  gl.useProgram(programInfo.program);

  gl.uniformMatrix4fv(
      programInfo.uniformLocations.projectionMatrix,
      false,
      projectionMatrix);
  gl.uniformMatrix4fv(
      programInfo.uniformLocations.modelViewMatrix,
      false,
      modelViewMatrix);

  gl.viewport(0, 0, gl.canvas.clientWidth, gl.canvas.clientHeight);
  gl.drawArrays(gl.POINTS, 0, 150);

  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.clearColor(1.0, 1.0, 1.0, 1.0);
  gl.clearDepth(1.0);
  gl.disable(gl.DEPTH_TEST);
  gl.disable(gl.BLEND);
  //gl.depthFunc(gl.LEQUAL);
  gl.clear(gl.COLOR_BUFFER_BIT);

  gl.bindBuffer(gl.ARRAY_BUFFER, buffers.quadvbo);
  gl.vertexAttribPointer(
        programInfo1.attribLocations.vertexPosition,
        2,
        gl.FLOAT,
        false,
        4 * 4,
        0);
  gl.enableVertexAttribArray(programInfo1.attribLocations.vertexPosition);

  gl.vertexAttribPointer(
        programInfo1.attribLocations.textureCoord,
        2,
        gl.FLOAT,
        false,
        4 * 4,
        2 * 4);
  gl.enableVertexAttribArray(programInfo1.attribLocations.textureCoord);

  gl.useProgram(programInfo1.program);
  gl.activeTexture(gl.TEXTURE0);
  gl.uniform1i(programInfo1.uniformLocations.textureSampler, 0);
  gl.bindTexture(gl.TEXTURE_2D, buffers.texture);

  gl.viewport(0, 0, gl.canvas.clientWidth, gl.canvas.clientHeight);
  gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

}

function initShaderProgram(gl, vsSource, fsSource) {
  const vertexShader = loadShader(gl, gl.VERTEX_SHADER, vsSource);
  const fragmentShader = loadShader(gl, gl.FRAGMENT_SHADER, fsSource);

  const shaderProgram = gl.createProgram();
  gl.attachShader(shaderProgram, vertexShader);
  gl.attachShader(shaderProgram, fragmentShader);
  gl.linkProgram(shaderProgram);

  if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
    alert('Unable to initialize the shader program: ' + gl.getProgramInfoLog(shaderProgram));
    return null;
  }

  return shaderProgram;
}

function loadShader(gl, type, source) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader, source);
  gl.compileShader(shader);

  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    alert('An error occurred compiling the shaders: ' + gl.getShaderInfoLog(shader));
    gl.deleteShader(shader);
    return null;
  }

  return shader;
}
