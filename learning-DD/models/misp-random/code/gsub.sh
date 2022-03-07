for file in src/dd/* src/learning/* include/dd/* include/learning/*
do
sed -i -e 's/IndepSetInst/IndepSetInst2/g' $file
done
