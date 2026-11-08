'use strict';
importScripts('statistics.js');
self.onmessage = function(event) {
  const {id, values, n, k}=event.data;
  try { self.postMessage({id,result:MatrixStatistics.scan(new Float64Array(values),n,k)}); }
  catch(error) { self.postMessage({id,error:error.message}); }
};
