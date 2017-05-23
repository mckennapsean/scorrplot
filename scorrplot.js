// global variables for testing
var dataset
var vectors = []
var standardizedVectors = []

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
      vector[i] = parseFloat(dataset[row]['X' + (i + 1)])
    }
    vectors[row] = vector
  }
  console.log('vectors computed')

  // standardize vectors (mean-centered, scaled)
  for (let i = 0; i < vectors.length; i++) {
    // initialize standardized vector
    standardizedVectors[i] = []
    // compute mean
    let mean = 0
    for (let j = 0; j < vectors[0].length; j++) {
      mean += vectors[i][j]
    }
    mean = mean / vectors[0].length
    // subtract mean
    for (let j = 0; j < vectors[0].length; j++) {
      standardizedVectors[i][j] = vectors[i][j] - mean
    }
    // compute length
    let length = 0
    for (let j = 0; j < standardizedVectors[0].length; j++) {
      length += standardizedVectors[i][j] * standardizedVectors[i][j]
    }
    length = Math.sqrt(length)
    // scale length to 1
    for (let j = 0; j < standardizedVectors[0].length; j++) {
      standardizedVectors[i][j] = standardizedVectors[i][j] / length
    }
    length = 0
  }
  console.log('data standardized')

  // select primary & secondary points of interest
  let primaryVector = 0
  let secondaryVector = 1

  // compute initial projection

  // select new primary point

  // interpolate projections
})
