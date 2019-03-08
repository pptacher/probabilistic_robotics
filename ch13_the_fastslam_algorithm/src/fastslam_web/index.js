var sio = require('socket.io');
var ss = require('socket.io-stream');
path = require('path');
const fs = require('fs');
var forwarded = require('forwarded-for');
var debug = require('debug')('fastslam:io');

process.title = 'fastslam-io';

var port = process.env.PORT || 3000;

module.exports.listen = function(app) {
  io = sio.listen(app);

  debug('listening on *:' + port);

  io.total = 0;
  io.on('connection', function(socket){
    var fastslam = require('./fastslam');
    var tid=-1, src;
    var req = socket.request;
    var ip = forwarded(req);
    debug('New connection from ' + ip.ip);
    debug('total: '+(++io.total));
    //var ip = req.connection.remoteAddress;
    debug('client ip %s',ip.ip);
    debug('client ip %s', req.connection.remoteAddress);

    ss(socket).on('setup-stream', function(stream) {

      debug('setup-stream');
      [tid, src] = fastslam.instream();
      src.pipe(stream);
    })

    socket.on('disconnect', function(){
      if (tid!=-1) {
        fastslam.killthread(tid);
        tid=-1;
      }
      debug('user disconnected.');
      --io.total;

    });

  }); //io.on('connection', function(socket)

  return io;
}
