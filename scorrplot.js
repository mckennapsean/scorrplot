// global variables for testing
var dataset
var vectors = []
var standardizedVectors = []
var projection = []

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

  // calculate secondary projection vector
  let secondaryVectorProjection = []
  // dot product of primary & secondary
  let dotProduct = 0
  for (let j = 0; j < vectors[0].length; j++) {
    dotProduct += standardizedVectors[secondaryVector][j] * standardizedVectors[primaryVector][j]
  }
  // compute projection vector (y)
  let secondaryVectorProjectionLength = 0
  for (let j = 0; j < vectors[0].length; j++) {
    secondaryVectorProjection[j] = standardizedVectors[secondaryVector][j] - dotProduct * standardizedVectors[primaryVector][j]
    secondaryVectorProjectionLength += secondaryVectorProjection[j] * secondaryVectorProjection[j]
  }
  secondaryVectorProjectionLength = Math.sqrt(secondaryVectorProjectionLength)
  // standardize projection vector
  for (let j = 0; j < vectors[0].length; j++) {
    secondaryVectorProjection[j] /= secondaryVectorProjectionLength
  }
  console.log('calculated new plane')

  // compute initial projection
  // for each vector, insert into projection the [x,y] coordinate as list
  for (let i = 0; i < vectors.length; i++) {
    let currentProjection = [0, 0]
    // compute the dot-product of each vector for x and linear combination for y
    for (let j = 0; j < vectors[0].length; j++) {
      currentProjection[0] += standardizedVectors[i][j] * standardizedVectors[primaryVector][j]
      currentProjection[1] += standardizedVectors[i][j] * secondaryVectorProjection[j]
    }
    // update primary vector to be aligned exactly on the far right
    if (i === primaryVector) {
      currentProjection[0] = 1
      currentProjection[1] = 0
    // update secondary vector to map to outer edge of projected circle
    } else if (i === secondaryVector) {
      let currentX = currentProjection[0]
      let newY = Math.sqrt(1 - currentX * currentX)
      currentProjection[1] = newY
    }
    // save the vector projection
    projection[i] = currentProjection
  }
  console.log('projected vectors')

  // initialize svg and plot elements
  // 200 x 200 grid to easily support scaling of strokes
  let plot = d3.select('.content').append('svg').attr('viewBox', '-100 -100 200 200')
  let scorrplot = plot.append('g').attr('class', 'scorrplot')
  scorrplot.append('circle').attr('r', 100).attr('class', 'border')
  scorrplot.append('g').attr('class', 'points')

  // visualize projection

  // select new primary point

  // interpolate projections
})
