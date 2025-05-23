
set -e
set -o pipefail

# Use the OpenPBTA bucket as the default.
URL=${SPLICE_URL:-https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/pbta-splicing}
RELEASE=${SPLICE_RELEASE:-v11}

# Remove old symlinks in data
find data -type l -delete

# The md5sum file provides our single point of truth for which files are in a release.
curl -k --create-dirs $URL/$RELEASE/md5sum.txt -o data/$RELEASE/md5sum.txt -z data/$RELEASE/md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < data/$RELEASE/md5sum.txt | cut -d ' ' -f 2` release-notes.md)


# Download the items in FILES if not already present
for file in "${FILES[@]}"
do
  if [ ! -e "data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl -k $URL/$RELEASE/$file -o data/$RELEASE/$file
  fi
done


GENCODE39="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz"
cd data
if [ ! -e ${GENCODE39##*/} ]
then
  echo "Downloading ${GENCODE39##*/}"
  curl -k -O $GENCODE39
fi

# Check the md5s for everything we downloaded except CHANGELOG.md
cd $RELEASE
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt
cd ../../


# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done
