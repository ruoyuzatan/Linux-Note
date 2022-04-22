#1.安装FRP服务（云服务器和本地都需要安装）
https://github.com/fatedier/frp/releases #下载链接
tar -zxvf viroblast.tar.gz #解压
#2.云服务器端配置
cd frp_0.23.0_linux_amd64
#赋予frps可执行权限
chmod +x frps
#将frps放到path目录下
cp frps /usr/bin/
  
#建立配置文件目录
mkdir /etc/frp/
#将配置文件复制到对应的配置文件目录
cp frps.ini /etc/frp/

#开机启动
vim /etc/systemd/system/frps.service
#文件内容
[Unit]
Description=FRP Server Daemon

[Service]
Type=simple
ExecStartPre=-/usr/sbin/setcap cap_net_bind_service=+ep /usr/bin/frps
ExecStart=/usr/bin/frps -c /etc/frp/frps.ini
Restart=always
RestartSec=20s
User=nobody
PermissionsStartOnly=true

[Install]
WantedBy=multi-user.target
#保存

#FRPS配置

sudo vim frps.ini

[common]
bind_port = 6001 #服务器端与本地通讯的端口
bind_udp_port = 7003 #点对点通讯用的端口
token = zhiwushengli #密钥

max_pool_count=5 #最大连接数


vhost_http_port = 6080 #客户端映射的端口，即通过该端口访问http网页

dashboard_port = 6100 #控制面板端口
dashboard_user = ruoyu
dashboard_pwd = 9536581999
#保存

#服务器端口配置
#去服务器端口管理页面，自定义，开放tcp端口：6000/6100

#3.本地配置
cd frp_0.23.0_linux_arm64
#赋予frps可执行权限
chmod +x frpc
#将frpc放到path目录下
cp frpc /usr/bin/

#建立配置文件目录
mkdir /etc/frp/
#将配置文件复制到对应的配置文件目录
cp frpc.ini /etc/frp/

#FRPC配置（客户端1，被访问的主机）
sudo vim frpc.ini

[common]
server_addr = 39.107.101.102 #云服务器地址
server_port = 6001 #云服务与本地通讯端口
token = zhiwushengli #密钥

[ssh]
type = tcp #协议
local_ip = 202.198.133.159 #本地内部ip
local_port = 22 #本地ssh端口，见ssh设置
remote_port = 6022 #远程控制端口，即在ssh应用里填的端口

[p2p_ssh]
type = xtcp
#只有sk一致的用户才能访问到此服务
sk = hanlei
local_ip = 202.198.133.159
local_port = 22

[httpname]
type = http
local_port = 80 #本地web服务端口
local_ip = 127.0.0.1 #如果custom_domains用外网ip，则注释掉
custom_domains = 域名(需解析外网ip）/外网ip

[plugin_socks] #将服务器作为网页代理服务器
type = tcp
remote_port = 6066
plugin = socks5
plugin_user = ruoyu
plugin_passwd = 9536581999
user_encryption = true
user_compression = true

#自启动设置

vim /etc/systemd/system/frpc.service

#文件内容
[Unit]
Description=FRP Client Daemon
After=network.target
Wants=network.target

[Service]
Type=simple
ExecStart=/usr/bin/frpc -c /etc/frp/frpc.ini
Restart=always
RestartSec=20s
User=nobody

[Install]
WantedBy=multi-user.target
#保存

#FRPC配置（客户端2，访问的主机）
[common]
server_addr = 39.107.101.102
server_port = 6001
token = zhiwushengli

[p2p_ssh_visitor]
type = xtcp
# xtcp 的访问者
role = visitor
# 要访问的 xtcp 代理的名字
server_name = p2p_ssh
sk = hanlei
# 绑定本地端口用于访问 ssh 服务
bind_addr = 192.168.31.164
bind_port = 7004

#重新启动服务
systemctl start frps #启动
systemctl stop frps #停止
systemctl restart frps #重启
systemctl status frps #查看状态
systemctl enable frps #开机启动frp

systemctl start frpc
systemctl stop frpc
systemctl restart frpc
systemctl status frpc
systemctl enable frpc

