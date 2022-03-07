for file in ./dd/*
do
sed -i -e 's/Edge2Node2MsgPass/Edge2NodeMsgPass/g' $file
done
for file in ./learning/*
do
sed -i -e 's/Edge2Node2MsgPass/Edge2NodeMsgPass/g' $file
done
