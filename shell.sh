## 登录
git config --global user.name "silhouette99"
git config --global user.email "13321218978@163.com"
## 关联仓库密钥
ssh-keygen -t ed25519 -C "13321218978@163.com"
cat -al /home/pzz/.ssh/id_ed25519_git.pub
## 测试链接git
ssh -T git@github.com
## 添加密钥
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519_git
## 
git clone git@github.com:silhouette99/Analysis_of_data_on_glioma_and_neuron_interactions.git


git init
git remote add origin git@github.com:silhouette99/Analysis_of_data_on_glioma_and_neuron_interactions.git
git pull origin main
git add .
git commit -m "first push"
git push origin main


git init
git remote add origin git@github.com:silhouette99/Analysis_of_data_on_glioma_and_neuron_interactions.git
git pull origin main
git add .
git commit -m "second push"
git push origin main

