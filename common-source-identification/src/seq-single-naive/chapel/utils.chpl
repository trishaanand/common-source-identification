/* Use sorting. */
use Sort;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Add a directory to a file name. */
proc addDirectory(fileName : string, dir : string) : string {
  return dir + "/" + fileName;
}

/* Get an array of the file names of the images in directory dir. */
proc getImageFileNames(dir : string) {
    var imageFiles = listdir(dir);
    sort(imageFiles);
    return addDirectory(imageFiles, dir);
}

/* Write a real array to a file. */
proc write2DRealArray(array : [] real, fileName :string) {
  assert(array.rank == 2);
  var file = open(fileName, iomode.cw);
  var channel = file.writer();

  for i in array.domain.dim(1) {
    for j in array.domain.dim(2) {
      if (i > j) {
        channel.writef("%{#####.#####} ", array[i, j]);
      } else if (i < j) {
        channel.writef("%{#####.#####} ", array[j, i]);
      } else {
        channel.writef("%{#####.#####} ", 0.0);
      }
    }
    channel.writeln();
  }
}