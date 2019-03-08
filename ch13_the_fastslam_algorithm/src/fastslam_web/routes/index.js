var express = require('express');
var router = express.Router();

/* GET home page. */
router.get('/', function(req, res, next) {
  //res.render('index', { title: 'Fastslam' });
  res.render('index0');
});


module.exports = router;
