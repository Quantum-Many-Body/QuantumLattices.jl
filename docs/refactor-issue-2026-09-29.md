# 重构 Issue：Frameworks 层概念重构（FrameworkElement / 删 Frontend / 删 Action / @delegate / 契约显式化）

日期：2026-09-29
范围：QuantumLattices（主战场）、TightBindingApproximation、ExactDiagonalization、QuantumClusterTheories（同步迁移）
性质：**破坏性变更**，各包同步发 breaking 版本。

## 背景与动机

现行 Frameworks 层（`src/Frameworks.jl`）的类型层级被用作"能力打包器"而非语义分类：`Assignment`/`Algorithm` 继承 `LatticeModel` 导致本体错乱（计算节点 is-a 物理模型）和类型边界洞（TBA 的 `H<:LatticeModel` 可接受 Assignment/Algorithm）；`Frontend` 是不声明任何接口的空抽象层；`Action` 的协议（`==`/`update!`）经查证为从未被下游实现的死代码（全部 11 个具体 Action 均无 `update!`，`Assignment.map` 仅服务于该死通道）；`valtype` 的语义沿继承链漂移，且 Assignment/Algorithm 上的定义经全仓库核查无任何消费者。本次重构只做**概念层**修正：不搬移文件、不修已知 bug（见文末"不在范围"）。

## 目标类型树

```
FrameworkElement                     # 新顶层：参数协议 + 命名/持久化/show + valtype 派生机制
├── LatticeModel                     # 物理模型（语义层；valtype 义务在此层）
│   ├── Formula / Generator(Static/Categorized/Operator)
│   └── （下游：TBA / ED / ImpuritySolver / SCMF 等"前端角色"实现）
├── Algorithm{L<:LatticeModel, P, M} # 执行上下文：model + parameters + map + timer
└── Assignment{A, P, T, D}           # 计算节点：action + parameters + dependencies + data
                                     # （A 无界；map 字段删除，7→6 字段）

Action：类型删除，角色保留（规格 struct 无父类，靠 run!/datatype dispatch）
Data：保留（Tuple(data) 协议 + datatype 的 <:Data 断言）
```

## 决议明细

### A1. 引入 `FrameworkElement` 顶层抽象

- 新增 `abstract type FrameworkElement end`；`LatticeModel <: FrameworkElement`；`Algorithm`、`Assignment` 直接 `<: FrameworkElement`（不再是 LatticeModel）。
- 以下方法的签名从 `::LatticeModel` 上移到 `::FrameworkElement`（只改签名不改实现）：
  - `Parameters` 的 `@generated` fallback（现 `:409-415`）
  - `contenttoconfig`（`:418-425`）/ `contenttocache`（`:454-461`）
  - `show` 默认实现（`:376-378` 区域）
  - `config`/`stamp`/`basename`/`pathof`/`str`（`:432-516`）
  - `qlsave(target::Symbol, ...)`（`:518`）
  - `valtype` 实例级转发（`:388`）、`scalartype`/`eltype` 通用派生（`:396-406`）——**机制**上移；"子类必须实现 valtype"的**义务**声明留在 LatticeModel 层 docstring。
- 通用 `==`/`isequal`（`:376` 区域）留在 LatticeModel 层（域收缩为真正的模型）。
- `Action`/`Data` 不入树。不新增通用 `update!` 默认实现。
- docstring 草案：`FrameworkElement` = "Frameworks 子模块管理的全部元素的抽象父类：晶格模型、算法、赋值。提供统一协议：参数（Parameters/update!）、命名（str/basename/pathof）、持久化（config/stamp/qlsave）、显示，以及可选的 valtype 协议（实现类型级 valtype 则自动获得 scalartype/eltype）"。

### A2. 删除 `Frontend` 与 `Action`（硬删，无别名过渡）

- QuantumLattices：
  - 删 `abstract type Frontend`（`:1114`）与 export（`src/QuantumLattices.jl:50`）。
  - `Algorithm{F<:Frontend,...}` → `Algorithm{L<:LatticeModel,...}`（`:1255`）；`datatype` 的 F 边界同步（`:1292`）。
  - **字段改名**：`Algorithm.frontend` → `Algorithm.model`（类型参数 F→L）。下游所有 `.frontend` 访问机械替换。
  - 删 `abstract type Action`（`:1121`）、其 `==`/`isequal`（`:1122-1123`）、`update!` no-op（`:1130`）。
  - `Assignment{A}`：A 无界；**删除 `map` 字段**（`:1159`）；构造器去 map 参数（`:1164`）；注册入口去 map 位置参数（`:1375-1388`）；`update!(assign)` 简化为只更新自身 `parameters`（`:1176-1182`）：
    ```julia
    function update!(assign::Assignment; parameters...)
        length(parameters)>0 && (assign.parameters = update(assign.parameters; parameters...))
        return assign
    end
    ```
  - `datatype`/options 机器（`:1192-1239`）的 `A<:Action` 边界全部解除；callable 中 `action::Action` 改为无界。
  - test/Frameworks.jl 的示例类型同步改造（`TBA{F} <: LatticeModel`、action struct 无父类、注册调用去 map）。
- 下游父类替换：`TBA`（TBA/Core.jl:274）、`ED`（ED/Core.jl:267）、`ImpuritySolver`（QCT/Core.jl:49）改为 `<: LatticeModel`。
- TBA 的 H 边界**不动**（`H<:LatticeModel`；嵌套前端是活概念，实例为 MeanFieldTheory 的 `SCMF <: TBA{K, N<:PureTBA{K}, Nothing}`）。`matrix`/`dimension` 的嵌套分支签名 `<:Frontend` → **`<:TBA`**（TBA/Core.jl:314、338）。
- HatsugaiKohmotoModels 为 stale 包（引用已不存在的 AbstractTBA/TargetSpace），不在本次范围。
- 不加 trait 函数；"前端/action"作为角色写进文档。

### A3. valtype 语义 + 接口透明性原则 + `@delegate` 宏

- 删 `valtype(::Type{<:Assignment})`（`:1168`，无消费者）；Assignment 无 valtype 概念。
- 保留 `valtype(::Type{<:Algorithm{L}}) = valtype(L)`（改指向 `model`）；`scalartype`/`eltype` 经 FrameworkElement 层通用派生自动获得；下游 ED 的 `scalartype(::Type{<:Algorithm{<:ED}})`（ED/Core.jl:280）可删。
- **接口透明性原则**写入 Algorithm docstring：Algorithm 是其 model 的 transparent proxy，`f(alg,...) ≡ f(alg.model,...)`，仅注入执行上下文（timer）。
- 实现并导出 `@delegate` 宏（Algorithm 专用），规格：
  ```julia
  @delegate function matrix(tba::TBA, k=nothing; gauge=:icoordinate) ... end   # 实例级
  @delegate inject=(:timer,) function eigen(tba::TBA, k; timer=tbatimer, o...) ... end  # 注入
  @delegate function kind(::Type{<:TBA{K}}) where K ... end                     # 类型级
  @delegate @inline f(x::X) = body                                              # 短形式+宏组合
  ```
  - 展开 = 原方法 + 转发方法。实例级：第一参数 `x::X`（X<:LatticeModel）→ `f(x::Algorithm{<:X}, ...) = f(x.model, ...)`；类型级：`f(::Type{<:X}, ...)` → `f(::Type{<:Algorithm{L}}, ...) where {L<:X} = f(L, ...)`。
  - **实施修正（2026-09-29）**：类型级转发的 `where` 必须带 `L<:X` 界（spec 原文为无界 `where L`）。无界时不同下游包对同名函数（如 TBA 与 ED 各自的 `kind`）生成完全相同的转发签名，多包同载互相 method overwrite；加界后签名不相交。
  - 转发方法沿用原参数名、默认值、定义形式（长/短）、包装宏（未知宏一律复制到转发版；`@generated` 不支持，docstring 注明）。
  - `inject=(...)`：同名 Symbol 约定；声明的 kwarg 在转发签名中保留但默认值改为 `x.同名字段`（默认覆盖语义），并显式传入转发调用；宏期校验 Symbol 必须是原签名真实 kwarg；仅实例级。
  - 解析管线：剥 `where` 层 → 剥包装 macrocall → 识别 `function`/`=` 本体。
  - 边界（docstring 写明）：只标注纯查询接口；`update!`/`Parameters` 永不标注；`@delegate` 必须是最外层宏；docstring 放更外层。
- 下游迁移：TBA 的 `kind` Union 惯用法（Core.jl:286-288）、`dimension`/`matrix` 转发（:295,:321-324）、`eigen/eigvals/eigvecs` 三组样板（:351-398）、ED 的 prepare!/release!/matrix/eigen 转发（Core.jl:280-360）改用 `@delegate`。
- 卫生项：TBA 内部 `datatype(::Type{D}, ::Union{Nothing,AbstractVector})`（Core.jl:341-342）改名（如 `matrixeltype`），消除与 QL 导出 `datatype` 的同名歧义。

### B1. `datatype`：显式注册可选，推断兜底保留

- 保持二元签名 `datatype(::Type{A}, ::Type{F}) where {A, F<:LatticeModel}`（A 无界）。
- 显式注册 = 下游定义更具体方法（dispatch 优先），**可选**；现存 11 个 Action 零迁移。
- 兜底推断（`Core.Compiler.return_type`）保留；`@assert` 报错改为指导性文字："datatype 推断失败：请显式定义 `datatype(::Type{X}, ::Type{Y}) = YourData`，或检查对应 `run!` 的类型稳定性"。

### B2. `dependencytypes`：依赖形状声明

- 新增导出函数：`dependencytypes(::Type) = nothing`（默认不校验）；下游 opt-in，如 `dependencytypes(::Type{<:GroundStateExpectation}) = (EDEigen,)`。
- 只支持定长；校验在注册入口汇流处（`:1381`），构造期执行：个数 + 逐位 `dep isa Assignment{<:T}`；报错写明期望与实得。
- 重构后（非本次）：已声明 Action 的 run! 内 `@assert isa(dependencies,...)` 冗余可清理。

### B3. 定义性 options 上移 Action 字段

- 迁移清单：
  - TBA `DensityOfStates`：收 `fwhm/ne/emin/emax`；`InelasticNeutronScatteringSpectra`：收 `fwhm/rescale`（`check` 留 hints）。
  - ED `EDEigen`：收 `nev/which`。
  - QCT `DynamicalSpectra`：收 `η/rescale`。
  - 构造器改为关键字字段形式；run! 改为读 `assign.action` 字段而非 options；`options(...)` 声明同步收缩。
- 留在 options（方法性 hints）：`gauge/infinitesimal/tol/maxiter/krylovdim/v₀/verbosity/release/showinfo/check`。
- 白名单机器（`options/hasoption/checkoptions/optionsinfo`）保留。
- **options 契约**写入文档：options = 执行提示，不参与 Assignment 身份，不触发缓存重算；需要不同定义的输出 → 构造新 Assignment。

### C1. 对偶 callable 保留 + 文档澄清

- `(alg)(assign)` 与 `(assign)(alg)` 均保留为公开 API。
- docstring 写明："**括号内的对象是参数权威**"；依赖递归走反向的原因（依赖继承调用点上下文）一并写明。
- 微清理：`checkoptions::Bool` 从位置参数改为关键字参数（`:1343`、`:1355`、递归处 `:1349/:1361`）。

### C2. 缓存契约成文 + 构造期参数键校验

- **缓存契约**（写入文档）：每个 Assignment 是单槽缓存；命中条件 = `isdefined(:data)` 且存储参数与当前上下文参数在 `(atol, rtol)` 内容差匹配；`action`/`dependencies` 为 const 不参与比较；options 不触发重算；不做递归依赖身份比较（`X.data` 是自洽快照，命中时与依赖当前状态无关；未命中时各节点局部检查保证正确）。
- 推论写明：① 参数容差内视为相同；② 强制重算的唯一途径是构造新 Assignment；③ 共享依赖在不同参数上下文间切换会抖动重算，需并存则复制依赖实例。
- 构造期校验：注册入口（`:1382`）本地参数键必须 ⊆ `alg.parameters` 键，未知键报错（修掉惰性多余键陷阱与反向 `match` 的 KeyError 边角）。
- `match` 保留 `Base.match` 扩展，不改名。

### C3. 字段分组正式化 + 等值规则

| 组 | Algorithm | Assignment | 消费者 |
|---|---|---|---|
| 标识 | `name` | `name` | `==`、`show`、`str` |
| 存储 | `dir` | `dir` | 仅 pathof/持久化 |
| 计算 | `model`、`parameters`、`map` | `action`、`parameters`、`dependencies` | `==`、`show`、config/stamp |
| 观测 | `timer` | — | 仅计时 |
| 结果 | — | `data` | 读取/落盘/绘图/`==`（带保护） |

- `Algorithm` 的 `==`/`isequal`：比较 `(name, model, parameters, map)`——**dir 移出**（行为变化：同模型不同目录判等）。
- 新增 `Assignment` 的 `==`/`isequal`：比较 `(name, action, parameters, dependencies)` + `data`（带 undef 保护）：
  ```julia
  d₁, d₂ = isdefined(a₁, :data), isdefined(a₂, :data)
  d₁ == d₂ || return false
  !d₁ && return true
  return a₁.data == a₂.data
  ```
  action 字段比较显式走 `efficientoperations`（Action 删除后无类型可分派）；独立 action 相等比较由下游按需一行定义。
- `isequal` 同构。

## 实施顺序与验收

1. **Phase 1 — QuantumLattices**：A1 → A2 → A3（含 @delegate）→ C3 → B1/B2 → C1/C2 微清理；`test/Frameworks.jl` 同步改造；`Pkg.test()` 全绿。
2. **Phase 2 — TightBindingApproximation**：父类/字段名/`@delegate` 迁移 + B3 options 上移；测试全绿。
3. **Phase 3 — ExactDiagonalization**：同上（EDEigen 收 nev/which）；测试全绿。
4. **Phase 4 — QuantumClusterTheories**：同上（DynamicalSpectra 收 η/rescale）；测试全绿。
5. 四包联调。

每 Phase 完成标准：该包 `Pkg.test()` 通过，无遗留对已删除符号（`Frontend`/`Action`/`frontend` 字段/`Assignment.map`）的引用。

## 不在本次范围（重构完成后另行处理）

**已知 bug（推迟修复）**：
- TBA `matrix` 嵌套分支丢 k（TBA/Core.jl:338-340）——影响 MeanFieldTheory SCMF 活路径，优先级最高；
- TBA `FermiSurface` 权重计算漏传 options（Core.jl:929）；
- ED v₀ 路径死代码/坏签名（ED/Core.jl:159、195、199）；
- ED `SectorFilter` 置零不剔除（Core.jl:258）；GreenFunction 体系未纳入 Assignment 框架。

**不搬移文件/依赖外迁**：HDF5/TimerOutputs/Latexify 外迁、dlmsave/持久化移 extension、Boundary/Embedding 归置——均不在本次范围。

---

# 附：bug 修复决议（2026-09-29 讨论定稿）

- **Bug 1（TBA matrix 嵌套分支丢 k，Core.jl:341）**：修复——转发 k（`matrix(getcontent(tba,:H), k; gauge=gauge, infinitesimal=infinitesimal)`），gauge/infinitesimal 维持外层默认下传。
- **Bug 2（TBA FermiSurface 权重漏传 options，Core.jl:935）**：修复——`eigvecs(tba, k; options...)`。全包扫描确认唯一漏网点。
- **Bug 3（ED v₀ 路径，Core.jl:159/:188-199）**：修复——① 向量 v₀ 走 KrylovKit 4 位置参数形式（Int v₀ 路径不变）；② 多 sector 版归一化后统一用 `m.ket`（Sector 键）查找，未知键报错；docstring 对齐。
- **Bug 4（ED SectorFilter 置零，Core.jl:258）**：**最终结论：代码维持元素级投影原状不动**，仅 docstring 写明设计契约——"元素级线性投影 + `OperatorSum` 经 `add!` 保证不含零元素（公共 API 不变量）= 集合级有效剔除"。讨论时的"幽灵零本征值"担忧经实证不复现（`add!` 丢零，QuantumOperators.jl:597-598）；曾尝试的集合级剔除方法经两轮复审认定为冗余（行为与不变量保证的结果完全一致），已撤除。净产出：docstring 契约 + 一组钉住该不变量链的回归测试（过滤结果无零块、元素级投影行为、正能量谱 sector 选择）。
- **Bug 5（ED GreenFunction 未纳入 Assignment 框架）**：**不动，关闭**。结论：GF 是交互式求值对象，与 Assignment 的"一次性计算+缓存结果"是两种合法设计，不强行套框架；RawStderrLogger 保留。

---

# 附：文档讨论产生的追加改名决议（2026-09-29 晚）

- **Assignment 字段 `action` → `task`**：action 角色改名为 task（computation task）。理由：物理包中 "action" 与作用量撞词；且改名后 `Assignment` 获得字面意义——assign a task to an algorithm。无类型牵连（Action 类型已删，角色词无 Base.Task 混淆风险，prose 首次出现写全 computation task）。波及：QL + 下游各包 run! 内的 `.action` 访问、docstring、文档。
- **Algorithm 字段 `model` → `frontend`（撤回 A2 决议②）**：理由：字段应按角色命名（与 task 对称），所装对象是"体系在具体算法下的预备形态"（TBA/ED 实例），叫 model 抹掉物理模型与算法预备形态的区别；frontend 的编译器语义准确，且角色词在文档中有定义，不再是化石。波及：全部迁移时 `.frontend`→`.model` 的改动机械还原。
- **Algorithm 存在理由的文档表述**（教程 7.1 采用）：执行层横向统一（存储/恢复/计时）+ 高一级参数管理（map：物理参数→terms 参数）+ Assignment 工厂与缓存权威。

---

# 附：文档更新决议（2026-09-29 逐节讨论定稿）

## 总原则

- 章节结构保留；"frontend/task"作为**角色词**（无类型、小写、无 @ref 链接），Data/Algorithm/Assignment/FrameworkElement 保留类型与链接。
- **FrameworkElement 在教程中完全不点名**（含 7.8），只在 man/Frameworks.md 出现；教程统一用 "share the model interface" 表述。
- 教程不教显式 datatype 注册（只讲推断路径）；教 dependencytypes 与 @delegate。
- 数值旋钮的缺省占位一律用 **NaN** 而非 nothing（类型稳定），推广到所有涉及包与测试。
- 所有 `[Frontend](@ref)`/`[Action](@ref)` 断链清扫。

## ToyTBA.jl 终版要点

`ToyTBA{M<:LatticeModel,T} <: LatticeModel`，字段 `hamiltonian`+`table`；写全契约四件套 valtype/Parameters/update!/contenttoconfig + show；`matrix` 两方法标 `@delegate`；EigenSystem/DensityOfStates 无父类；DensityOfStates 收字段 `emin::Float64=NaN, emax::Float64=NaN, ne::Int=101, σ::Float64=0.1`（关键字构造器），run! 用 `isnan` 判断并从 `assignment.task` 解构；`dependencytypes(::Type{<:DensityOfStates}) = (EigenSystem,)`；EigenSystem 保留 showinfo option（文中点明其为纯执行提示）。

## 第 7 章逐节决议

- **7.1**：引言段与 minimal-example 段不动；五角色压缩为三对象（model/Algorithm/Assignment）+两附属（task/Data）；动机改为从具体算法着手（TBA 要单粒子二次型、ED 要占据数表象稀疏矩阵 → 每算法一个 LatticeModel 子类做适配，同时统一 chap6 的表示）；末尾段不提 FrameworkElement，只说 "share the model interface"；流程图 `Frontend(model)`→`ToyTBA(model)`。
- **7.2**：标题不动；ToyTBA 代码块换新形态（含 valtype 行）；matrix 块标 @delegate 并新增接口透明性段落（transparent proxy，仅注入执行上下文）。
- **7.3**：Haldane 示例全段不动；"algorithm is a LatticeModel" 表述两处改 "shares the model interface"；:229 段改为"Algorithm 单一具体类型 vs frontend 角色每包各填"；结尾段 action→task、"assignment is a task"→"assigns a computation task to the algorithm"。
- **7.4**：标题改 "Task, Data and `run!`"；EigenSystem 去父类 + "无父类因 dispatch 于具体 struct"教学；run! 中 `.action`→`.task`、`matrix(algorithm.frontend, k)`→`matrix(algorithm, k)`；datatype 段只讲推断。
- **7.5**：对偶调用段加 C1 规则"参数权威是括号内的对象"；DOS 代码块按 ToyTBA 终版重写；删除 Assignment map 的括号句；注册段加"per-call 参数键须 ⊆ 算法参数键"。
- **7.6**：开头新增定义性 vs 提示的分类段（缓存契约为理由）；DOS options 声明删除；typo 示例改用 `shwoinfo`；optionsinfo 叙事改为"showinfo 来自依赖链"；checkoptions 改关键字参数表述；删 map 从句；结尾示例拆为 `update!`+重跑 与 `DensityOfStates(emin=-4.0, emax=4.0)` 新 assignment 两个演示。
- **7.7**：开篇句改 "share the model interface"；缓存小节补两句（options 不参与身份比较；强制重算=构造新 assignment）；共享依赖抖动不写进教程；其余不动。
- **7.8**：标题改 "Why the interface has this shape"；Julia 限制段保留但结论翻转；新增"为什么 frontend/task 没有抽象类型"段（标记抽象不强制任何东西；契约由 datatype/dependencytypes/类型化字段强制；继承是打包不是分类学）；表格按三对象+两角色重写（Assignment 无 map）；Data 论述保留、"nothing to inherit"措辞软化；FrameworkElement 不点名。
- **7.9/Summary**：表格不动；"reuses its actions"→"tasks"；Summary 逐条按新词汇改写（algorithm/assignment "share the model interface"）。

## 小处

- `1-introduction.md:17` "generic frontend"→"generic platform"。
- `6-latticemodel.md` 通读核查类型归属断言。
- `man/Frameworks.md` 引言重写：FrameworkElement 为根（此处点名）、角色清单更新、清单加 dependencytypes/@delegate。
- showinfo 保留并点明"纯执行提示"身份。

## 验收

Documenter docs build 必须通过（@example 块真实执行 ToyTBA）。实施顺序：先代码改名（action→task、model→frontend 回退，QL+下游五包测试全绿），再文档改写，最后 docs build。
