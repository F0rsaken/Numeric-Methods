var typedArray = new Uint8Array([0.1, 0.2, 0.7])

var trace1 = {
    type: 'bar',
    width: typedArray,
    x: [1, 2, 3],
    xaxis: 'x',
    y: [3, 2, 1]
}

var data = [trace1];

var layout = {
  xaxis: {'domain': [0, 1]},
};

Plotly.newPlot('div1', data, layout);