import os
import sys
import shutil

def replace_text_in_file(file_path, old_text, new_text, keep_original=False):
    """
    指定されたファイル内のテキストを置換する
    ファイル名に置換前のテキストが含まれている場合は、新しいファイルを作成して元のファイルは残す
    
    Args:
        file_path: 対象ファイルのパス
        old_text: 置換前のテキスト（例：runA）
        new_text: 置換後のテキスト（例：runB）
        keep_original: 元のファイルを残すかどうか（デフォルト：False）
    """
    try:
        # ファイルが存在するか確認
        if not os.path.isfile(file_path):
            print(f"エラー: ファイルが見つかりません: {file_path}")
            return False
        
        # ディレクトリとファイル名に分割
        directory = os.path.dirname(file_path)
        filename = os.path.basename(file_path)
        
        # ファイル名に置換対象のテキストが含まれているか確認
        if old_text in filename:
            # 新しいファイル名を作成
            new_filename = filename.replace(old_text, new_text)
            new_file_path = os.path.join(directory, new_filename)
            
            # ファイルをコピー（元のファイルは残す）
            shutil.copy2(file_path, new_file_path)
            print(f"ファイルをコピー: {file_path} → {new_file_path}")
            
            # コピーしたファイルの内容を置換
            with open(new_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            new_content = content.replace(old_text, new_text)
            
            with open(new_file_path, 'w', encoding='utf-8') as f:
                f.write(new_content)
            
            print(f"成功: {new_file_path}")
            print(f"ファイル内容の '{old_text}' を '{new_text}' に置換しました")
            print(f"元のファイル '{file_path}' は残されています")
            return True
        else:
            # ファイル名に置換対象が含まれていない場合は、そのまま内容を置換
            if keep_original:
                # 元のファイルを残したい場合
                base, ext = os.path.splitext(filename)
                new_filename = f"{base}_{new_text}{ext}"
                new_file_path = os.path.join(directory, new_filename)
                shutil.copy2(file_path, new_file_path)
                target_file = new_file_path
                print(f"ファイルをコピー: {file_path} → {new_file_path}")
            else:
                target_file = file_path
            
            with open(target_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            new_content = content.replace(old_text, new_text)
            
            if content == new_content:
                print(f"警告: '{old_text}' はファイル内に見つかりませんでした")
                return False
            
            with open(target_file, 'w', encoding='utf-8') as f:
                f.write(new_content)
            
            print(f"成功: {target_file}")
            print(f"'{old_text}' を '{new_text}' に置換しました")
            return True
        
        return True
    except Exception as e:
        print(f"エラー: {e}")
        return False


if __name__ == "__main__":
    # コマンドライン引数から取得する場合
    if len(sys.argv) == 4:
        file_path = sys.argv[1]
        old_text = sys.argv[2]
        new_text = sys.argv[3]
        replace_text_in_file(file_path, old_text, new_text)
    else:
        # 直接指定する場合はここを編集してください
        print("使用方法:")
        print("  python changeRunNumber.py <ファイルパス> <置換前> <置換後>")
        print()
        print("例:")
        print("  python changeRunNumber.py config.txt runA runB")
        print()
        print("機能:")
        print("  - ファイル名に置換対象が含まれている場合は、新しいファイルを作成")
        print("  - 元のファイルは残される")
