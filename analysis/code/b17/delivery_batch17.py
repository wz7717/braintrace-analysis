"""Create the B17 delivery inside the controller's existing clock and byte ledger.

This module does not run a comparator or read features, models, caches, or B13
results. Importing it has no packaging side effect.
"""
from pathlib import Path
import hashlib
import json
import zipfile

BLOCK = 64 * 1024
PLAN_SHA256 = 'af571e775d42164bee031749856c4f6377f9f4792a9d383491190e11d0dac315'


class _DeliveryStopped(Exception):
    pass


def package_delivery(run_root, preparation_root, delivery_root, receipt_path,
                     state, ledger, global_reason):
    """Copy new evidence only; return stop_reason and preserve any partial output.

    Call after the run archive and before the controller's final budget record.
    The caller retains responsibility for its authoritative terminal run files.
    """
    run_root, preparation_root = Path(run_root), Path(preparation_root)
    delivery_root, receipt_path = Path(delivery_root), Path(receipt_path)
    archive_path = delivery_root.with_name(delivery_root.name + '.zip')
    records, stop_reason, charged_bytes = {}, None, 0

    def check():
        nonlocal stop_reason
        stop_reason = stop_reason or global_reason()
        if stop_reason:
            raise _DeliveryStopped(stop_reason)

    def charge(path, count, current_size):
        nonlocal charged_bytes
        ledger.note_write(path, count, 'output')
        # note_write charges the delta; sample() must see the cumulative size.
        ledger.current[Path(path)] = current_size
        charged_bytes += count

    def source_chunks(path):
        check()
        with path.open('rb') as source:
            while True:
                check()
                chunk = source.read(BLOCK)
                if not chunk:
                    break
                yield chunk

    def write_chunks(path, chunks):
        check()
        path.parent.mkdir(parents=True, exist_ok=True)
        digest, size = hashlib.sha256(), 0
        with path.open('wb') as target:
            for chunk in chunks:
                for offset in range(0, len(chunk), BLOCK):
                    check()
                    piece = chunk[offset:offset + BLOCK]
                    count = target.write(piece)
                    target.flush()
                    size += count
                    charge(path, count, size)
                    digest.update(piece[:count])
                    check()
        return {'bytes': size, 'sha256': digest.hexdigest()}

    def put(name, data):
        record = write_chunks(delivery_root / name, (data,))
        records[name] = {'file': name, **record}
        return records[name]

    def put_json(name, value):
        check()
        return put(name, json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False).encode('utf-8'))

    def copy(path, name):
        record = write_chunks(delivery_root / name, source_chunks(path))
        records[name] = {'file': name, **record}

    class ArchiveWriter:
        """Append-only ZIP sink: every actual write is charged, including headers."""
        def __init__(self, target):
            self.target, self.size = target, 0
            self.digest = hashlib.sha256()

        def tell(self):
            return self.size

        def seek(self, *args):
            raise OSError('The budgeted ZIP sink is append-only')

        def flush(self):
            self.target.flush()

        def write(self, data):
            check()
            for offset in range(0, len(data), BLOCK):
                check()
                piece = data[offset:offset + BLOCK]
                count = self.target.write(piece)
                self.target.flush()
                self.size += count
                charge(archive_path, count, self.size)
                self.digest.update(piece[:count])
                check()
            return len(data)

    try:
        check()
        delivery_root.mkdir(parents=True, exist_ok=False)
        summary_path = run_root / 'results/SUMMARY.json'
        summary = json.loads(b''.join(source_chunks(summary_path)).decode('utf-8-sig')) if summary_path.exists() else {}
        is_real = state.get('mode') == 'real'
        kind = summary.get('input_data_kind', 'REAL_GSE25219_FROZEN_FINAL_PLM' if is_real else 'SYNTHETIC_B17_END_TO_END')
        put_json('TERMINAL_LOCK_20261108.json', {
            'schema': 'braintrace.batch17.terminal_execution_lock.v1',
            'input_data_kind': kind, 'authorized': is_real, 'enabled': False,
            'real_execution_authorized': is_real, 'real_data_execution_enabled': False,
            'no_automatic_continuation_retry_or_follow_on': True,
            'run_state_at_package_start': state.get('status'),
            'global_stop_at_package_start': state.get('global_stop'),
            'method_stops_at_package_start': state.get('method_stops', {}),
            'authoritative_terminal_run_state': str(run_root / 'RUN_STATE.json'),
            'authoritative_final_budget': str(run_root / 'FINAL_BUDGET.json'),
            'packaged_state_is_snapshot': True})
        for path in sorted(Path(__file__).parent.iterdir()):
            check()
            if path.is_file() and path.suffix.lower() in {'.py', '.r', '.md'}:
                copy(path, 'code/' + path.name)
        # Preparation contains newly written metadata only. Referenced numerical
        # paths are copied as text; they are never followed or opened here.
        for path in sorted(preparation_root.rglob('*')):
            check()
            if path.is_file() and '__pycache__' not in path.parts and path.suffix.lower() in {'.json', '.csv', '.md', '.txt'}:
                copy(path, 'preparation/' + path.relative_to(preparation_root).as_posix())
        for name in ('OBSERVATIONS.json', 'ROSTER.csv'):
            check()
            if (run_root / name).exists():
                copy(run_root / name, 'run/' + name)
        put_json('run/RUN_STATE_SNAPSHOT.json', state)
        for path in sorted((run_root / 'results').rglob('*')):
            check()
            if path.is_file() and path.suffix.lower() in {'.json', '.csv'}:
                copy(path, 'run/results/' + path.relative_to(run_root / 'results').as_posix())
        put_json('PRIVATE_ARTIFACT_LOCATORS_20261108.json', {
            'authoritative_run_root': str(run_root),
            'cache_reference_matrices_and_RDS_not_copied': str(run_root / 'cache'),
            'query_cache_not_copied': str(run_root / 'queries'),
            'job_receipts_and_full_logs_private': str(run_root / 'jobs'),
            'run_bundle_not_duplicated': str(run_root / 'RUN_BUNDLE.zip'),
            'authoritative_final_budget_written_after_delivery': str(run_root / 'FINAL_BUDGET.json'),
            'original_failures_NOT_RUN_and_partial_outputs_retained': True})
        put_json('VERSION_BRIDGE_20261108.json', {
            'schema': 'braintrace.batch17.version_bridge.v1',
            'file_path_basis': 'Relative old entry paths resolve from this delivery directory',
            'B15_frozen_plan_sha256_reused': PLAN_SHA256,
            'B15_plan_and_protocol_entry': '../BrainTrace_R8_Batch15_External_Comparator_Freeze_20261108/INDEX.html',
            'B16_synthetic_implementation_entry': '../BrainTrace_R8_Batch16_External_Adapter_20261108/INDEX.html',
            'B13_historical_results_entry_only': '../BrainTrace_R8_Batch13_External_Evaluation_20261108/INDEX.html',
            'B13_RAW_donor_equal_Top1_historical_percent_reference_only': 89.44,
            'old_B13_result_files_read_or_recomputed': False,
            'old_assets_models_environments_copied_or_rehashed': False,
            'no_new_statistics_tuning_training_or_advantage_claim': True})
        methods = summary.get('methods', {})
        lines = []
        for method, result in methods.items():
            scope = result.get('scopes', {}).get('ALL', {})
            values = scope.get('contributing_donor_equal', {})
            lines.append(f"| {method} | {result.get('completion_status')} | {values.get('top1')} | {values.get('top3')} | {values.get('mrr')} |")
        table = '\n'.join(lines) or '现有SUMMARY不存在；不生成方法点估计。'
        readme = f'''# 第十七批 冻结final PLM外部比较

观察来源：`{kind}`。本包复制本次run已生成的原精度九类观察、原名单和现成导出；下表直接摘取SUMMARY，不重新计算方法分数或端点。None表示现成导出未给点估计。

| 方法 | 完成状态 | donor-equal Top1 | Top3 | MRR |
| --- | --- | --- | --- | --- |
{table}

真实运行使用B15接受的final PLM来源绑定，参数与端点继承冻结计划；B16仅为合成实现验证。本次属于已见B13结果之后的描述性比较，不追溯称为预注册，不作CI/P值、重采样、调参、重训或优势判断。九类参与排名不等于九区域外部验证；准确率范围与原N/D以现成SUMMARY和ROSTER为准。单方法失败/NOT_RUN及不利结果原样保存，不以成功子集替代全队列。

B13 RAW donor-equal Top1 89.44%仅作历史索引引用；未读取或重算B13数值文件。旧B15/B16/B13定位在VERSION_BRIDGE中，不改写旧包或身份。没有复制真实特征、cache矩阵、大RDS、模型或旧环境；完整job日志、参考/cache及原run ZIP保留在PRIVATE_ARTIFACT_LOCATORS记录的私有位置。

TERMINAL_LOCK的enabled=false，表示本次授权结束后不会自动续跑。preparation中的准入锁和run状态是历史snapshot，不能替代当前终止锁。原run `{run_root / 'RUN_STATE.json'}` 与稍后收口的 `{run_root / 'FINAL_BUDGET.json'}` 权威；本包和ZIP不在最终预算写入后回写。打包复制、写入与ZIP增长仍在同一总clock/ledger内逐64KiB监控、删除不退款；wall/RSS是采样soft monitor，不冒称OS硬限或零overshoot。已有全规模资源与失败状态只按本次run记录说明。

本包为本地私有交付，不授予数据或方法再分发权，不自动发布、push或联系第三方。未修改稿件、旧科学结论、环境或冻结B15方案。
'''
        put('README_20261108.md', readme.encode('utf-8'))
        targets = ['README_20261108.md', 'TERMINAL_LOCK_20261108.json',
                   'run/OBSERVATIONS.json', 'run/ROSTER.csv', 'run/RUN_STATE_SNAPSHOT.json',
                   'PRIVATE_ARTIFACT_LOCATORS_20261108.json', 'VERSION_BRIDGE_20261108.json', 'DELIVERY_MANIFEST.json']
        if summary_path.exists():
            targets.append('run/results/SUMMARY.json')
        links = ''.join(f'<li><a href="{name}">{name}</a></li>' for name in targets)
        page = ('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>BrainTrace Batch17</title>'
                '<h1>第十七批 冻结外部比较</h1><p>本次授权执行的交付snapshot；真实执行开关已终止，原run最终状态与预算权威。</p><ul>' + links + '</ul></html>')
        put('INDEX.html', page.encode('utf-8'))
        put_json('DELIVERY_MANIFEST.json', {'batch': 17, 'input_data_kind': kind,
                 'file_path_basis': 'Relative files resolve from this delivery directory',
                 'files': [records[name] for name in sorted(records)],
                 'hashes_cover_new_delivery_files_only': True,
                 'state_and_budget_authoritative_at_original_run': str(run_root),
                 'packaged_state_is_snapshot': True})
        check()
        missing = [name for name in targets if not (delivery_root / name).exists()]
        if missing:
            raise FileNotFoundError(missing)
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with archive_path.open('wb') as target:
            sink = ArchiveWriter(target)
            # The non-seekable sink makes ZIP headers/descriptors append-only.
            # On cap stop, closing attempts hit check() before further writes.
            with zipfile.ZipFile(sink, 'w', zipfile.ZIP_DEFLATED) as package:
                for name in sorted(records):
                    check()
                    with package.open(delivery_root.name + '/' + name, 'w') as member:
                        for chunk in source_chunks(delivery_root / name):
                            check()
                            member.write(chunk)
            check()
        receipt = {'batch': 17, 'input_data_kind': kind,
            'status': 'LOCAL_DELIVERY_SNAPSHOT_TERMINAL_REAL_DISABLED',
            'entry': str(delivery_root / 'INDEX.html'), 'delivery_root': str(delivery_root),
            'manifest': {**records['DELIVERY_MANIFEST.json'], 'file': str(delivery_root / 'DELIVERY_MANIFEST.json')},
            'zip': {'file': str(archive_path), 'bytes': sink.size, 'sha256': sink.digest.hexdigest()},
            'method_completion_statuses_from_existing_summary': {name: value.get('completion_status') for name, value in methods.items()},
            'state_at_package_start': state.get('status'), 'method_stops': state.get('method_stops', {}),
            'packaged_state_is_snapshot': True, 'current_new_entry_targets_exist': True,
            'authoritative_terminal_run_state': str(run_root / 'RUN_STATE.json'),
            'authoritative_final_budget_written_after_delivery': str(run_root / 'FINAL_BUDGET.json'),
            'package_bytes_charged_before_receipt_write': charged_bytes,
            'real_execution_authorized': is_real, 'real_data_execution_enabled': False,
            'no_automatic_continuation': True}
        # The receipt is outside the ZIP to avoid a cyclic ZIP/receipt identity.
        write_chunks(receipt_path, (json.dumps(receipt, ensure_ascii=False, indent=2,
                                              allow_nan=False).encode('utf-8'),))
        check()
        return {'stop_reason': None, 'entry': str(delivery_root / 'INDEX.html'),
                'receipt': str(receipt_path), 'archive': str(archive_path),
                'bytes_charged': charged_bytes, 'packaged_state_is_snapshot': True}
    except _DeliveryStopped:
        return {'stop_reason': stop_reason, 'partial_delivery_root': str(delivery_root),
                'partial_archive': str(archive_path), 'bytes_charged': charged_bytes,
                'regular_packaging_stopped': True, 'partial_files_preserved': True}
