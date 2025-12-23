import time

def print_timestamp(message):
    """
    在消息前添加时间戳，并打印到终端。

    参数:
    - message: 要打印的消息。
    """
    # 获取当前时间戳
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    
    # 构建带时间戳的消息
    full_message = f"[{timestamp}] {message}"
    
    # 打印消息到终端
    print(full_message)
