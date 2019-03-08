const binding = require('bindings')('addon');
const { Readable } = require('stream');
var stream = require('stream');

module.exports.instream = function () {
  const inStream = new Readable({
    objectMode: true,
    read(chunk) {}
  });
  //inStream._read = function(chunk) {};
  var tid = binding.startThread(inStream.push.bind(inStream));
  return [tid, inStream];

}

module.exports.killthread = function (tid) {

  binding.killThread(tid);
  return;

}
