// global variables for testing
var dataset
var vectors = []
var svectors = []

d3.csv('data/test.csv', (e, d) => {
  console.log('data loaded')

  // store data in browser
  dataset = d

  // figure out how many dimensions of data
  let dimensions = Object.keys(dataset[0]).length - 2

  // compute vectors and store them
  for (let row in dataset) {
    let vector = []
    for (let i = 0; i < dimensions; i++) {
      vector[i] = dataset[row]['X' + (i + 1)]
    }
    vectors[row] = vector
  }
  console.log('vectors computed')

  // standardize vectors (mean-centered, scaled)
  for (let i = 0; i < vectors.length; i++) {
    // compute mean
    let mean = 0
    for (let j = 0; j < vectors[0].length; j++) {
      mean += vectors[i][j]
    }
    mean = mean / vectors[0].length
    // subtract mean
    for (let j = 0; j < vectors[0].length; j++) {
      svectors[i][j] = vector[i][j] - mean;
    }
    // scale length to 1
    for (let j = 0; j < vectors[0].length; j++) {}
  }
  console.log('data standardized')

  // select primary & secondary points of interest
  let ppoi = 0
  let spoi = 1
})
