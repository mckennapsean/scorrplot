var dataset
var vectors = []

d3.csv('data/test.csv', (e, d) => {
  console.log('data loaded')
  // console.log(d)

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
})
